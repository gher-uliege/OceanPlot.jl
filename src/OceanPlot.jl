module OceanPlot

if VERSION >= v"0.7"
    using Statistics
else
    using DataArrays
end

using PyPlot
using PyCall
using PyCall: PyObject, pyimport
@static if VERSION < v"0.7"
    using MAT
end
using NCDatasets
using DIVAnd
using LinearAlgebra

if VERSION < v"0.7"

    @pyimport numpy.ma as ma
    @pyimport matplotlib.patches as pypatch


    masked(S,mask) = pycall(ma.array, Any, S, mask=mask)


    # Plotting DataArrays (masked arrays for python)

    pyma(S) = masked(S.data, S.na)
    PyPlot.pcolor(z::DataArray; kws...) = pcolor(pyma(z); kws...)
    PyPlot.pcolor(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
    pcol(z::DataArray; kws...) = pcolor(pyma(z); kws...)
    pcol(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
    PyPlot.plot(x::DataArray, y::DataArray; kws...) = plot(pyma(x), pyma(y); kws...)

    # Plotting using NaNs
    NaNpyma(S) = masked(S, isnan.(S))
    pcol(z::Array{T,2}; kws...) where T = pcolor(NaNpyma(z); kws...)
    pcol(x,y,z::Array{T,2}; kws...) where T = pcolor(x,y,NaNpyma(z); kws...)
else

    function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
        numpy_ma = PyCall.pyimport("numpy")["ma"]
        pycall(numpy_ma["array"], Any, coalesce.(a,zero(T)), mask=ismissing.(a))
    end

    PyObject(a::Adjoint) = PyObject(copy(a))
    PyObject(a::Transpose) = PyObject(copy(a))
end

export plot_coastline, pcol, listfiles, set_aspect_ratio, patch, plotvecstd

patch(x,y; kwargs...) = gca().add_patch(pypatch.Polygon(cat(2,x,y); kwargs...))


function ncview(fname,varname,slide)
    Dataset(fname) do ds
        pcol(ds[varname][slide...]'; cmap="jet")
    end
end



"""
Plots the coastline from the file `fname`. The file `fname` is a .mat file with the variables `ncst` and `Area`.
"""

function plot_coastline(fname = joinpath(ENV["HOME"],"Data","Coastline","gshhs_l.mat");
                         patchcolor = [.8,.8,.8], linewidth = 2.,zorder = nothing,plottype = :plot)
    if fname in ["l","i","c","h"]
        fname = joinpath(ENV["HOME"],"Data","Coastline","gshhs_$(fname).mat")
    end

    xl = xlim()
    yl = ylim()

    c = matread(fname)
    k = round.(Int64,c["k"])
    ncst = c["ncst"] :: Array{Float64,2}
    area = c["Area"][:,1] :: Vector{Float64}
    index = find(area .> 0);
    ax = gca();

    cplot(ncst) =
        if plottype == :plot
            plot(ncst[:,1],ncst[:,2],"k-",linewidth = linewidth);
        else
            ax.add_patch(pypatch.Polygon(ncst,color = patchcolor, zorder = zorder))
        end

    for l=1:length(index)
        i = index[l];
        j = k[i]+1:k[i+1]-1;
        cplot(ncst[j,:])
    end

    xlim(xl)
    ylim(yl)

    nothing
end


function plotmap(bathname = joinpath(ENV["HOME"],"projects","Julia","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_4.nc"),
                  patchcolor = [.8,.8,.8], coastlinecolor = nothing)

    xl = xlim()
    yl = ylim()
    # work-around
    xl = xl[1]:0.1:xl[2]
    yl = yl[1]:0.1:yl[2]

    bx,by,b = DIVAnd.extract_bath(bathname,true,xl,yl)
    if patchcolor !== nothing
        contourf(bx,by,b', levels = [-1e5,0],colors = [patchcolor])
    end

    if coastlinecolor !== nothing
        contour(bx,by,b', levels = [-1e5,0],colors = coastlinecolor, linestyles = "-")
    end
end


function hview(fname, varname, subindex...; orientation="horizontal", cmap="jet",
               lonname = "lon", latname = "lat", vmin = nothing, vmax = nothing)
    ds = Dataset(fname);
    S = nomissing(ds[varname][subindex...],NaN);
    lon = nomissing(Dataset(fname)[lonname][subindex[1]],NaN);
    lat = nomissing(Dataset(fname)[latname][subindex[2]],NaN);
    close(ds)

    pcolor(lon,lat,copy(S'); cmap=cmap, vmin = vmin, vmax = vmax)
    set_aspect_ratio()
    colorbar(orientation=orientation)
    plotmap()
    xlim(extrema(lon)...)
    ylim(extrema(lat)...)
end

"""
Fixes the aspect ratio of a plot.
"""
function set_aspect_ratio()
    ax = gca()
    as = cos(mean([ylim()...]) * pi/180)
    ax.set_aspect(1/as)
end



function plotvecstd1(x,y,u1,v1; scale = 1, scaleu = scale,
                     scalev = scale, scalestd = 1)
    um = mean(u1)
    vm = mean(v1)

    up = u1 - mean(u1)
    vp = v1 - mean(v1)

    P = Symmetric([up⋅up  up⋅vp ; 0  vp⋅vp ], :U) / (length(up)-1)

    λ, U = eig(P)
    xy = [U[i,1] * sqrt(λ[1]) * cos(θ) + U[i,2] * sqrt(λ[2]) * sin(θ) for θ = linspace(0,2π,100), i = 1:2];

    patch(x + scaleu*(scalestd*xy[:,1] + um),y + scalev*(scalestd*xy[:,2] + vm);
          color = "lightblue", zorder = 1)

    quiver([x],[y],[scaleu*um],[scalev*vm]; angles = "xy", scale_units = "xy", scale = 1, zorder = 2)

end

function plotvecstd(x,y,u,v; scale = 1, scaleu = scale,
                    scalev = scale, scalestd = 1, mincount = 10,
                    legendpos = [], legendvec = [1,0], legendcolor = "r")

    for j = 1:size(x,2)
        for i = 1:size(x,1)
            u1 = u[i,j,:]
            v1 = v[i,j,:]
            if sum(.!ismissing.(u1)) >= mincount && sum(.!ismissing.(v1)) >= mincount
                plotvecstd1(x[i,j],y[i,j],u1.data,v1.data;
                            scaleu = scaleu, scalev = scalev,
                            scalestd = scalestd)
            end
        end
    end

    if length(legendpos) == 2
        quiver(legendpos[1],legendpos[2],scaleu*legendvec[1],scalev*legendvec[2];
               angles = "xy", scale_units = "xy", scale = 1, zorder = 2,
               color = legendcolor)
    end

end


"""
List all files starting from `topdir` with the provided `extension`.
"""
function listfiles(topdir = "."; extension = "")
    list = String[]

    for (root,dirs,files) in walkdir(topdir)
        for file in files
            if length(extension) == 0
                push!(list, joinpath(root, file))
            else
                if endswith(file,extension)
                    push!(list, joinpath(root, file))
                end
            end
        end
    end
    return list
end

end # module
