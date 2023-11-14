function [im,nx,ny] = dicom_thumb(info)

imx  = info.IconImageSequence.Item_1.PixelData;
im   = typecast(imx,'uint8');
nx   = 64;
ny   = 64;
im   = reshape(im,nx,ny);