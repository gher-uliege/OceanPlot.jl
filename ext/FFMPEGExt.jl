module FFMPEGExt

using FFMPEG
using Printf

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

end
