module OceanPlot

using NullableArrays
using DataArrays
using PyPlot
using PyCall
using MAT

@pyimport numpy.ma as ma

export plot_coastline, pcol, listfiles, set_aspect_ratio

# Plotting NullableArrays (masked arrays for python)

na2pyma(S) =  pycall(ma.array, Any, S.values, mask=S.isnull)
PyPlot.pcolor(x,y,z::NullableArray; kws...) = pcolor(x,y,na2pyma(z); kws...)
PyPlot.pcolor(z::NullableArray; kws...) = pcolor(na2pyma(z); kws...)
pcol(x,y,z::NullableArray; kws...) = pcolor(x,y,na2pyma(z); kws...)
pcol(z::NullableArray; kws...) = pcolor(na2pyma(z); kws...)
PyPlot.plot(x::NullableArray, y::NullableArray; kws...) = plot(na2pyma(x), na2pyma(y); kws...)

# Plotting DataArrays (masked arrays for python)

pyma(S) =  pycall(ma.array, Any, S.data, mask=S.na)
PyPlot.pcolor(z::DataArray; kws...) = pcolor(pyma(z); kws...)
PyPlot.pcolor(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
pcol(z::DataArray; kws...) = pcolor(pyma(z); kws...)
pcol(x,y,z::DataArray; kws...) = pcolor(x,y,pyma(z); kws...)
PyPlot.plot(x::DataArray, y::DataArray; kws...) = plot(pyma(x), pyma(y); kws...)

pcol{T}(z::Array{T,2}; kws...) = pcolor(pycall(ma.array, Any, z, mask=isnan.(z)); kws...)
pcol{T}(x,y,z::Array{T,2}; kws...) = pcolor(x,y,pycall(ma.array, Any, z, mask=isnan.(z)); kws...)



"""
Plots the coastline from the file `fname`. The file `fname` is a .mat file with the variables `ncst` and `Area`.
"""
function plot_coastline(fname = joinpath(ENV["HOME"],"Data","Coastline","gshhs_l.mat"))
    xl = xlim()
    yl = ylim()

    c = matread(fname)
    k = round(Int64,c["k"])

    index = find(c["Area"] .> 0);

    for l=1:length(index)
        i = index[l];
        j = k[i]+1:k[i+1]-1;

        #  patch(c.ncst[j,1],c.ncst[j,2],[.8 .8 .8]);
        plot(c["ncst"][j,1],c["ncst"][j,2],"k-",linewidth=2.);
    end

    xlim(xl)
    ylim(yl)
end

"""
Fixes the aspect ratio of a plot.
"""
function set_aspect_ratio()
    ax = gca()
    as = cos(mean([ylim()...]) * pi/180)
    ax[:set_aspect](1/as)
end

"""
List all files starting from `topdir` with the provided `extension`.
"""
function listfiles(topdir = "."; extension = "")
    list = []

    for (root,dirs,files) in walkdir(".")
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
