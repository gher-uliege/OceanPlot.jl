module OceanPlot

using FFMPEG
using Statistics
using PyPlot
using PyCall
using PyCall: PyObject, pyimport
using MAT
using NCDatasets
using DIVAnd
using LinearAlgebra
using Printf

# allow for plotting with missing values
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

export plot_coastline, pcol, listfiles, set_aspect_ratio, patch, plotvecstd

const mpl = PyPlot.matplotlib

patch(x,y; kwargs...) = gca().add_patch(mpl.patches.Polygon(cat(x,y,dims=2); kwargs...))


function ncview(fname,varname,slide)
    Dataset(fname) do ds
        pcolormesh(ds[varname][slide...]')
    end
end

function romsview(fname::AbstractString,varname,slide; kwargs...)
    Dataset(fname) do ds
        v = ds[varname]
        lon = NCDatasets.coord(v,"longitude")[:]
        lat = NCDatasets.coord(v,"latitude")[:]
        v = ds[varname][slide...]
        Δlon = lon[2,2] - lon[1,1]
        Δlat = lat[2,2] - lat[1,1]

        pcolormesh(lon .- Δlon/2,lat .- Δlat/2,v; kwargs...)
        set_aspect_ratio()

        return lon,lat,v
    end
end


function romsview(fnames::Vector,varname,slide; kwargs...)
    for fname in fnames
        lon,lat,v = romsview(fname,varname,slide; kwargs...)
        plot([lon[1,1],lon[end,1],lon[end,end],lon[1,end],lon[1,1]],
             [lat[1,1],lat[end,1],lat[end,end],lat[1,end],lat[1,1]],"k")
    end
    colorbar()

    bathname = joinpath(ENV["HOME"],"projects","Julia","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_1.nc")
    @show bathname
    plotmap(bathname)
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
    index = findall(area .> 0);
    ax = gca();

    cplot(ncst) =
        if plottype == :plot
            plot(ncst[:,1],ncst[:,2],"k-",linewidth = linewidth);
        else
            ax.add_patch(mpl.patches.Polygon(ncst,color = patchcolor, zorder = zorder))
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


function plotmap(bathname = joinpath(ENV["HOME"],"projects","Julia","DIVAnd-example-data","Global","Bathymetry","gebco_30sec_4.nc");
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


function hview(fname, varname, subindex...; orientation="horizontal", cmap=nothing,
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
    as = cosd(mean(ylim()))
    ax.set_aspect(1/as)
end



function plotvecstd1(x,y,u1,v1; scale = 1, scaleu = scale,
                     scalev = scale, scalestd = 1)
    um = mean(u1)
    vm = mean(v1)

    up = u1 .- mean(u1)
    vp = v1 .- mean(v1)

    P = Symmetric([up⋅up  up⋅vp ; 0  vp⋅vp ], :U) / (length(up)-1)

    λ, U = eigen(P)
    xy = [U[i,1] * sqrt(λ[1]) * cos(θ) + U[i,2] * sqrt(λ[2]) * sin(θ) for θ = LinRange(0,2π,100), i = 1:2];

    patch(x .+ scaleu*(scalestd*xy[:,1] .+ um),
          y .+ scalev*(scalestd*xy[:,2] .+ vm);
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
                plotvecstd1(x[i,j],y[i,j],u1,v1;
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


"""

Example:

```julia
varname = "temp"; cl = (6.216519f0, 16.935743f0); outname = "SST.mp4"
OceanPlot.animation(fname,varname,cl,outname)
```
"""
function animation(fname,varname,cl,outname; figprefix = replace(outname,".mp4" => "-"), framerate  = 5, qtile = (0.05,0.95))

    get(v,n) = (ndims(v) == 4 ? v[:,:,end,n] : v[:,:,n] )

    Dataset(fname) do ds
        figure(figsize=(10,6));

        v = ds[varname]

        if cl == nothing
            q = [quantile(skipmissing(vec(get(v,k))), qtile) for k = 1:size(v)[end]];
            cl = (minimum(getindex.(q,1)), maximum(getindex.(q,2)))
            @show cl
        end

        lon = NCDatasets.coord(v,"longitude")[:]
        lat = NCDatasets.coord(v,"latitude")[:]

        for n = 1:size(v)[end]
            clf()
            v_slice = get(v,n)
            pcolormesh(lon,lat,nomissing(v_slice,NaN));
            OceanPlot.set_aspect_ratio();
            clim(cl)
            colorbar()
            title("SST $(ds["ocean_time"][n])")
            figname = figprefix * @sprintf("%05d.png",n)
            @info "saving $figname"
            savefig(figname)
        end
    end

    FFMPEG.exe("-y",
               "-r",string(framerate),
               "-i",figprefix * "%*.png",
               "-qmax","20",
               "-b:v","50k",
               outname)
end
end # module
