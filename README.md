Simple example:
```julia
using BoucheTrous, EasyFITS

filename = "/Users/ferreol/Data/RawData/SPHERE/BoucheTrou/reduced_SPHER.2019-11-02T00:16:29.340_2639773_IRDIS_OBJECT_7.996176s_16f_DB_K12_IMAGE,DUAL,CORONOGRAPHY.fits"

hdu = read(FitsHeader, filename);
img = read(FitsFile, filename; ext=1);
w   = read(FitsFile, filename; ext=2);
hdr = read(FitsHeader, filename);


badpix = w .==0;

#IRDIS left panel
leftbox = BoucheTrous.get_box(badpix[1:1024,:,:])
bouchetrous!(view(img,leftbox,:),view(badpix,leftbox,:);maxiter=1000)

#IRDIS right panel
rightbox = BoucheTrous.get_box(badpix[1025:end,:,:]) + (1025,0)
bouchetrous!(view(img,rightbox,:),view(badpix,rightbox,:);maxiter=1000)

## Detect and interpolate spurious bad pixel using median filter
using StatsBase,ImageFiltering
function redetectbadpix!(data, weights)
	med = mapwindow(median,data,(3,3,1))
	medw = mapwindow(median,weights,(3,3,1))
	P = (med .- data ).*(medw .+ 1e-6)
	md = mad(P)
	badpix = .!((-10*md) .< P .<  (10*md))
	data[badpix] .= med[badpix]
	weights[badpix] .= 0
end 
```