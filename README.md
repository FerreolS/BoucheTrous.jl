Simple example:
```julia
using BoucheTrous, EasyFITS

filename = "/Users/ferreol/Data/RawData/SPHERE/BoucheTrou/reduced_SPHER.2019-11-02T00:16:29.340_2639773_IRDIS_OBJECT_7.996176s_16f_DB_K12_IMAGE,DUAL,CORONOGRAPHY.fits"

hdu = read(FitsHeader, filename);
img = read(FitsFile, filename; ext=1);
w   = read(FitsFile, filename; ext=2);
hdr = read(FitsHeader, filename);


badpix = w .==0;
bouchetrous!(img,badpix)

#IRDIS left panel
leftbox = BoucheTrous.get_box(badpix[1:1024,:,:])
bouchetrous!(view(img,leftbox,:),view(badpix,leftbox,:);maxiter=1000)

#IRDIS right panel
rightbox = BoucheTrous.get_box(badpix[1025:end,:,:]) + (1025,0)
bouchetrous!(view(img,rightbox,:),view(badpix,rightbox,:);maxiter=1000)


```