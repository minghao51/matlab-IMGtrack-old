function [I rect1 rect2 rect3] = CroppingIMG(input,frame, path)
% This would crop at certain input, save the cropping detail at the path of
% the image

figure,imshow(input{round(mean(frame))},[])
[I rect1] = imcrop;
[I rect2] = imcrop;
[I rect3] = imcrop;

dlmwrite(fullfile(path, 'predata.txt'), ...
    [rect1; rect2; rect3; frame frame]);

