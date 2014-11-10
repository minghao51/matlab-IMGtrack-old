function [BestRow, BestCol] = MaxCorrelation2(ImageIn, template);
%based on http://www.youtube.com/watch?v=Q-OzmDen4HU
% the final BestRow and BestCal seems to be at the top left corner of the
% template.

%  imshow(ImageIn,[]);
[hI,wI] = size(ImageIn);

% cd(path)
% Template = imread('template.png');
[hT,wT] = size(template);

cc = normxcorr2(template, ImageIn);
[RowCc, Colc] = size(cc);

[max_cc, imax] = max(abs(cc(:)));
[ypeak, xpeak] = ind2sub(size(cc),imax(1));

BestRow = ypeak - (hT -1);
BestCol = xpeak - (wT -1);
%  h = imrect(gca, [BestCol BestRow (wT-1) (hT -1)]);