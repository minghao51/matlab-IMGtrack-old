function ImageOut = Gaussian_fn(ImageIn, hsize, sigma)
%Acquired from http://stackoverflow.com/questions/2773606/gaussian-filter-in-matlab/2773651#2773651
%# Create the gaussian filter with hsize = [5 5] and sigma = 2
% G = fspecial('gaussian',[5 5],2);
G = fspecial('gaussian',[hsize hsize],sigma);
%# Filter it
ImageOut = imfilter(ImageIn,G,'same');