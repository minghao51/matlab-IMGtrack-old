

% %% Clearing and some default settings
clear all
close all
%%

%Load All functions
 addpath('C:\Users\Minghao\Desktop\research project\GUI\functions')

%  %Path to the folders continues ImageJ processed
% tiff.fname_khamis1_Cor = ['C:\Users\Minghao\Desktop\research project' ...
%     '\Lung Datasets\MING_PATIENTS\KHAMIS1\haste_cor'];
% 
% % tiff.fname_khamis1_Sag = ['C:\Users\Minghao\Desktop\research project' ...
% %     '\Lung Datasets\MING_PATIENTS\KHAMIS1\haste_sag'];

% %Reading the data
% [tiff.IMGkhamis1Cor tiff.Corkhamis1Spacing] = ...
%     ReadingTiff([tiff.fname_khamis1_Cor, '\ResliceofHASTEcor.tif']);
% % [tiff.IMGkhamis1Sag tiff.Sagkhamis1Spacing] = ...
% %     ReadingTiff(tiff.fname_khamis1_Sag, '\ResliceofHASTEsag.tif');

%loading file name and paths
[tiff.Name,tiff.pName] = uigetfile('*', 'Load data');
if tiff.pName == 0, return; end

%Reading the data
[tiff.IMG tiff.PixelSpacing] = ...
    ReadingTiff([tiff.pName, tiff.Name]);

%% Narrowing the Data volume
% Being resliced projections,there are far too many frames, it would be best
%to narrow down the search area for tumour with cropped image and range of
%frames
% in the case of , it is manually determined that the tumour is 
% around 
% ie  khamis1_Cor = frame 82 - 107
%     khamis1_Sag = frame 92 - 114
    tiff.frame = [79 107];
    tiff.NumberF = tiff.frame(2)- tiff.frame(1);
% Searching for the initial crop coordinates and FeatureExtraction Coordinates
% Read for crop coordinate data and
% This search the Image folder for predata.txt, that is supposed to
% recorded the relevent cropping detaisl

tiff.predataName = 'predata.txt';
if exist(fullfile(tiff.pName, tiff.predataName), 'file')
    
    % Read the File, and output it into matrix
    M = dlmread(fullfile(tiff.pName, tiff.predataName));
    
    % Updating crop coordinates and method's option
    tiff.cropcor1 = M(1,:);
    tiff.cropcor2 = M(2,:);
    tiff.cropcor3 = M(3,:);
    
    try
        tif.frame = M(4,:);
        tiff.ExtractionCor = M(5,:);
    end
else
    % File does not exist.
    %Provide Cropping at approached location
    [temp tiff.cropcor1 tiff.cropcor2 tiff.cropcor3] = ...
        CroppingIMG(tiff.IMG, tiff.frame, tiff.pName);
end

% Crop now
for cnt = tiff.frame(1) : tiff.frame(2)
    tiff.I1{cnt - tiff.frame(1)+1} = imcrop(tiff.IMG{cnt},tiff.cropcor1);
    tiff.I2{cnt - tiff.frame(1)+1} = imcrop(tiff.IMG{cnt},tiff.cropcor2);  % times 2 because of the orignal half int
end

%% Gaussian filter
for cnt = 1 : tiff.NumberF
    tiff.GI1{cnt} = Gaussian_fn(tiff.I1{cnt}, 3,2) ;
    tiff.GI2{cnt} = Gaussian_fn(tiff.I2{cnt}, 3,2) ;
end

%% Segmentation of Tumour
tiff.OtsuBinV = 3;

for cnt = 1 : tiff.NumberF
    tiff.LOGI1{cnt} = otsu(tiff.GI1{cnt},tiff.OtsuBinV)*tiff.OtsuBinV;
end

clear IMG
for k = 1 : tiff.NumberF
    IMG(:,:,k) = tiff.LOGI1{k};
end
figure,imshow3D(IMG)

% Crop, would have to account for the crop pixel differences between
% the initial I1 and I3
%     tiff.LOGI1{cnt} = imcrop(tiff.LOGI1{cnt}, [tiff.cropcor3(1)-tiff.cropcor2(1) ....
%     tiff.cropcor3(2)-tiff.cropcor2(2) tiff.cropcor3(3) tiff.cropcor3(4)]);
for cnt = 1 : tiff.NumberF
    tiff.LOGI1{cnt}  = tiff.LOGI1{cnt}(15:25,:);
end
  %% Segmentation of Diaphragm
for cnt = 1 : tiff.NumberF
    tiff.LOGI2{cnt} = otsu(tiff.GI2{cnt},tiff.OtsuBinV)*tiff.OtsuBinV;
end
% figure,imshow(tiff.LOGI1{1},[])
% figure,imshow(tiff.LOGI2{1},[])

%% Features Extraction
%   Some parameters that enable GUI control, enable all here
tiff.Erode_DilateV = 0;
tiff.Bwareaopen = 0;
tiff.Erode_DilateV = 0;
tiff.Imfill = 0;
tiff.DynamicUpV = 0;

[ tiff ] = FeatureExtractTif( tiff.LOGI1,tiff.LOGI2, tiff );


