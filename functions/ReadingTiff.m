function [ImgOut PixelSpacing] = ReadingTiff(fname)
% This read in a multiImg Tiff, and output a ImgOut - Image of the Tiff at
% particular frame in the format of (:,:,k) and PixelSpacing, the pixel
% resolution of that tif, being indentical format to dicom's pixelspacing

%Finding the info for the tiff file
info = imfinfo(fname);
num_images = numel(info);
%From the info, infer the pixel spacing 
 PixelSpacing = [1/info(1,1).XResolution,1/info(1,1).YResolution];
for k = 1:num_images;
    temp = imread(fname, k, 'Info', info);
    ImgOut{k} = temp ;
end