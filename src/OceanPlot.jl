module OceanPlot

using NullableArrays
using DataArrays
using PyPlot
using PyCall
using MAT
using NCDatasets

@pyimport numpy.ma as ma
@pyimport matplotlib.patches as pypatch

export plot_coastline, pcol, listfiles, set_aspect_ratio, patch, plotvecstd

patch(x,y; kwargs...) = gca()[:add_patch](pypatch.Polygon(cat(2,x,y); kwargs...))

masked(S,mask) = pycall(ma.array, Any, S, mask=mask)
# Plotting NullableArrays (masked arrays for python)


na2pyma(S) = masked(S.values, S.isnull)
PyPlot.pcolor(x,y,z::NullableArray; kws...) = pcolor(x,y,na2pyma(z); kws...)
PyPlot.pcolor(z::NullableArray; kws...) = pcolor(na2pyma(z); kws...)
pcol(x,y,z::NullableArray; kws...) = pcolor(x,y,na2pyma(z); kws...)
pcol(z::NullableArray; kws...) = pcolor(na2pyma(z); kws...)
PyPlot.plot(x::NullableArray, y::NullableArray; kws...) = plot(na2pyma(x), na2pyma(y); kws...)

# Plotting DataArrays (masked arrays for python)

pyma(S) = masked(S.data, S.na)
PyPlot.pcolor(z::DataArray; kws...) = pcolor(pyma(z); kws...)
PyPlot.pcolor(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
pcol(z::DataArray; kws...) = pcolor(pyma(z); kws...)
pcol(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
PyPlot.plot(x::DataArray, y::DataArray; kws...) = plot(pyma(x), pyma(y); kws...)


# Plotting using NaNs
NaNpyma(S) = masked(S, isnan.(S))
pcol{T}(z::Array{T,2}; kws...) = pcolor(NaNpyma(z); kws...)
pcol{T}(x,y,z::Array{T,2}; kws...) = pcolor(x,y,NaNpyma(z); kws...)


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
            ax[:add_patch](pypatch.Polygon(ncst,color = patchcolor, zorder = zorder))
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

"""
Fixes the aspect ratio of a plot.
"""
function set_aspect_ratio()
    ax = gca()
    as = cos(mean([ylim()...]) * pi/180)
    ax[:set_aspect](1/as)
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
            if sum(.!isna.(u1)) >= mincount && sum(.!isna.(v1)) >= mincount
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
