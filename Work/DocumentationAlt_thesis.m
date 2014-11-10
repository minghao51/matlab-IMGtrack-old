%% Documents
%Test case for programs in GUI: loading all files as usual, but stop
% before segmenting, would be changes to test all different segmentation
% methods

%% Clearing and some default settings
clear all
close all

% docked figures are easier to manage
set(0,'DefaultFigureWindowStyle','docked')

% Matlab complains a lot in dock format
iptsetpref('ImshowInitialMagnification', 'fit')
warning('off','Images:imshow:magnificationMustBeFitForDockedFigure'); %Suppress Docked Fit Warning
%  set(0,'DefaultFigureWindowStyle','normal')

%Load All functions
 addpath('C:\Users\Minghao\Desktop\research project\GUI\functions')
% %  supress all figure!!! for publishing for now
% set(0,'DefaultFigureVisible','off')

handles.cropcor1 = [51.51,94.51,75.98,51.98];
handles.cropcor2 = [23.51,201.51,84.98,52.98];

%% Reading the Data
% Using matlab UI to navigate into the responding folder

[fName,pName] = uigetfile('*', 'Load data');
handles.pName = pName; % to pass this directory for other function
if pName == 0, return; end
%         dicomlist = dir(fullfile(pName,'Images','*.dcm'));
handles.dicomlist = dir(fullfile(pName, '*'));  %Generating a list of filename for reading in that directory, based on the 1st letter on the name
handles.dicomlist(~strncmp({handles.dicomlist.name}, fName(1), 1)) = []; % this yank out those files that isn't start with the same 1st letter as selected file
% handles.dicominfo = dicominfo(fullfile(pName,handles.dicomlist(1).name));   %reading dicominfo from the very 1st file
% Sorting the image, so that the order of time/frame are consistent
% chronologically
[handles.fname inx]= sort_nat({handles.dicomlist.name});
handles.dicomlist = handles.dicomlist(inx);

for cnt = 1 : numel(handles.dicomlist)
     handles.data{cnt} = im2double(dicomread(fullfile(pName,handles.dicomlist(cnt).name))); % directory reading of dicom files
     handles.info{cnt}  = dicominfo(fullfile(pName,handles.dicomlist(cnt).name)) ;
end

%Setting up ImgFolder for items
handles.ImageFolder = pName;                                   %passing the directory path to textbox/display

%% Searching for the initial crop coordinates and FeatureExtraction Coordinates
% Read for crop coordinate data and
% This search the Image folder for predata.txt, that is supposed to
% recorded the relevent cropping detaisl

handles.predataName = 'predata.txt';
if exist(fullfile(pName, handles.predataName), 'file')
    % File exists.  Do stuff....
%     % Display that notice the user there is pre-existing crop region
%     txtInfo = sprintf('There is pre-existing crop region \n cropping and thresholding will generate new predata.txt');
%     set(handles.txtbox, 'string', txtInfo);
    
    % Read the File, and output it into matrix
    M = dlmread(fullfile(pName, handles.predataName));
    
    % Updating crop coordinates and method's option
    handles.cropcor1 = M(1,:);
    handles.cropcor2 = M(2,:);
%     handles.MethodV = M(3,1);
    try
    handles.ExtractionCor = M(3,:);
%     catch exception
%         figure, imshow(handles.I1{cnt})
%         [impixel 
    end
else
    % File does not exist.
    % Display that notice the user there is not pre-existing crop region
    txtInfo = sprintf('There is no pre-existing crop region \n will generate predata.txt after cropping and thresholding \n handles.cropcor1');
    set(handles.txtbox, 'string', txtInfo);
    
end

% hence after crop and imagesegmentation, write it to the file, though this would
% require sharing of pName, fullFileName as handles
dlmwrite(fullfile(handles.pName, handles.predataName), ...
    [handles.cropcor1; handles.cropcor2])





%% For Original Image
%This section would crop the relevent section, while aso preserving the I1
%and I2 original cropping for comparison later on.

%Cropping and Saving as two cells
for cnt = 1 : numel(handles.dicomlist)
    handles.I1{cnt} = imcrop(handles.data{cnt},handles.cropcor1);
    handles.I2{cnt} = imcrop(handles.data{cnt},handles.cropcor2);  % times 2 because of the orignal half int
end

for cnt = 1 : numel(handles.dicomlist)
    oriI1{cnt} = handles.I1{cnt};
    oriI2{cnt} = handles.I2{cnt};
end

%Depiciting the Image
figure, imhist(handles.I1{1}), title('Original Image Histogram')
figure, imshow( handles.I1{1},[])
   
   %% Gaussian filter
% applying gaussian filter to OriginalLocal histogram and CLAHE processed frame
% just to check the influence 

    for cnt = 1 : numel(handles.dicomlist)
        handles.Gdata{cnt} = Gaussian_fn(handles.data{cnt}, 3,2) ;
        handles.GI1{cnt} = Gaussian_fn(handles.I1{cnt}, 3,2) ;
        handles.GI2{cnt} = Gaussian_fn(handles.I2{cnt}, 3,2) ;
    end  
            
%% Combining Histogram


%% Test of multiple frame 
% to compare histogram varies between frames
% obviously should be kept within the range of number of frames
% It seems the original histogram stack on top one another too much, hence,
% to compare them, I guess this test would be best to be use on Local
% Histogram, as it has the best spread, however, it is actually only known
% after comparing all of them. Filters such as Gaussian, median is thought 
% to have no influence here.

test = [1 5 13 20 22 25 30];

% Well, apparently, it seems like from the whole dicom image that the
% suppression is caused by blood in aorta arch

%% Computing the number of Outlying peaks on histogram


%% Segmentation of Tumour

LineWidth = 2;
%% Segmentation method 2: Local Otsu method
% I would expect the same result, but doing local should have the benifit
% it doing without the OtsuBinV control, at least, for most cases. In the
% 1st frame, it is observed that LO has better result, though whether that
% is conslusive remain to be seen

% In both case, the Gaussian filtered image is the best, Histogram
% Localisation and CLAHE doesn't seems to aid in segmentation much, though,
% perhaps looking at more frame, more data sets is require for any
% conclusion

handles.OtsuBinV = 2;

for cnt = 1 : numel(handles.dicomlist)

    handles.LOGI1{cnt} = otsu(handles.GI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;

end

% going through to depict imshow on test frames, 
for i = 1 : 7
    h=figure('Name','Local Otsu method on Gaussian data'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian data'});
    [B,L] = bwboundaries(handles.LOGI1{test(i)},'noholes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', LineWidth)
    end
    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
   
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\OtsuSegmentation\' handles.dicomlist(i).name '.png']) ;
end

%% Segmentation method 3: MCWS method
%Using Marker controlled watershed segmentation on I1, otsu for
% diaphragm I2. note watershed works like bwlabel, ie it label the
% different region with different values approprially
%% MCWS - Marker Controlled Water Segmentation method.
%Segmentation method to extract a segmentation of tumour with different
%frames

% Control for Background mask value,Foreground mask value and Foreground Objects connected pixel regions control
MBGM = 0; bgmV = 0;
MFGM = 0; fgmV= 0;
imregionalmaxV = 8;

for cnt = 1 : numel(handles.dicomlist)
%     handles.dataT{cnt} = handles.data{cnt};
    handles.MCWSI1{cnt} = MCWS1(handles.GI1{cnt}, bgmV, MBGM, fgmV, MFGM, imregionalmaxV );
end

[handles] = RegionMerging(handles);



MCWS1_PlotEveryVitals(handles.GI1{test(1)}, bgmV, MBGM, fgmV, MFGM, imregionalmaxV );

for i = 1: 7
    L =MCWS1(handles.GI1{test(i)}, bgmV, MBGM, fgmV, MFGM, imregionalmaxV );
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    h = figure, imshow(Lrgb,[]);
    title('Colored watershed label matrix (Lrgb)')
    
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\MCWSSeg\' handles.dicomlist(i).name '.png'])
end

for i = 1 : 7
    %     h=figure('Name','Marker Controlled Water Segmentation on Gaussian data'),...
    %         imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({'Marker Controlled Water Segmentation on Gaussian data'});
    %     [B,L] = bwboundaries(handles.MCWSI1{test(i)},'noholes');
    %     hold on
    %     for k = 1:length(B)
    %         boundary = B{k};
    %         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', LineWidth)
    %     end
    %     hold off
    %     set(gca,'position',[0 0 1 1],'units','normalized')
    
    h = figure('Name','Marker Controlled Water Segmentation on Gaussian data'),...
        imshow(handles.MCWSI1{test(i)},[],'InitialMagnification',100), title({'Marker Controlled Water Segmentation on Gaussian data'});
    
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\MCWSSeg\AfterReg\' handles.dicomlist(i).name '.png']) ;
end






%% Segmentation method 4:Template Matching method


% Loading template
template = double(imread([handles.ImageFolder 'template.png']));
template = template(:,:,1);

for cnt = 1 : numel(handles.dicomlist)
    [BestRow, BestCol] = MaxCorrelation2(handles.I1{cnt}, template);
    % inserting template at the correct location of the coronal plane
    % of picture
    [hT,wT] = size(template);
    handles.TMI1{cnt} = false(size(handles.I1{1}));
    handles.TMI1{cnt}(BestRow : BestRow + hT - 1, BestCol : BestCol + wT- 1 ) = template > 0;
end

% %Plot of Original Frame with template overlay on top of it's frame
% for i = 1: 3
%     figure('Name',' Original'),...
%         imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
%     
%     green = cat(3, zeros(size(handles.I1{1})), ones(size(handles.I1{1})), zeros(size(handles.I1{1})));
%     hold on
%     h = imshow(green);
%     % Use  TM1 map as the AlphaData for the solid green image overlay on
%     % Original data
%     set(h, 'AlphaData', handles.TMI1{test(i)}*0.3)
%     hold off
%     %     figure('Name','Template Matching method on Original'),...
%     %         imshow(handles.TMI1{test(i)},[],'InitialMagnification',100), title({'Template Matching method on Original'})
% end


    
for i = 1: 7
    h=figure('Name',' template matching segmentation on Gaussian data'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'});
    [B,L] = bwboundaries(handles.TMI1{test(i)},'noholes');
        hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', LineWidth)
    end
    hold off
    
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\TemplateSeg\' handles.dicomlist(i).name '.png'])
end

%% Segmentation method 5:Seed Growing Method
% This should potentially benefit from Histogram Equalisation thanks to the
% better spread of pixels. To be look at depth....

% Check for Exist for substructure, handles.ExtractionCor, 
% as normal exist() isn't enough.
doesVarExist = true ;% 1 = True
try
    handles.ExtractionCor;
catch
    doesVarExist = false ; % 0 = False
end

%Case for doesVarExist for manual selection or not
if doesVarExist == false
    % region seed growing selection of starting position
    figure, imshow(handles.I1{1},[]);
    [x,y,vals] = impixel;
else
    [x,y,vals] = impixel(handles.I1{1},...
        handles.ExtractionCor(1),handles.ExtractionCor(2));
end

hist(handles.I1{1},50);
prompt = {'Enter range fromo ROI'};
reg_maxdist = str2double(inputdlg(prompt))

%Need some way of automation
% if reg_maxdist ==0
%     reg_maxdist{cnt} =  graythresh(handles.I1{cnt})/3;
%     handles.I1{cnt}=regiongrowing(handles.I1{cnt},x,y,reg_maxdist{cnt});
% else
    
    for cnt = 1 : numel(handles.dicomlist)
        %     [nelements,xcenters] = hist(handles.I1{cnt});
        %     dydx=diff(nelements);
        %     level = graythresh(handles.I1{cnt});
        handles.dataT{cnt} = handles.data{cnt};
%         BackG = im2bw(handles.I1{cnt}, graythresh(handles.I1{cnt}));
%         handles.RSI1{cnt}(BackG==0) = 0;
        handles.RSI1{cnt}=regiongrowing(handles.I1{cnt},x,y,reg_maxdist);
    end
% end

for i = 1: 3
    figure('Name',' Original'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
    figure('Name','Seed Growing Method on Original'),...
        imshow(handles.RSI1{test(i)},[],'InitialMagnification',100), title({'Seed Growing Method on Original'})
end

    
for i = 1: 7
    h=figure('Name',' Seed Growing on Gaussian data'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'});
    [B,L] = bwboundaries(handles.RSI1{test(i)},'noholes');
        hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth',LineWidth )
    end
    hold off
    
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\SeedSeg\' handles.dicomlist(i).name '.png'])
end


%% Segmentation of Diaphragm
for cnt = 1 : numel(handles.dicomlist)
    handles.LOGI2{cnt} = otsu(handles.GI2{cnt},handles.OtsuBinV)*handles.OtsuBinV;
end
    

%% Features Extraction
% Planning to make a function for features extraction, verifying, plotting

%   Some parameters that enable GUI control, enable all here
handles.Erode_DilateV = 1;
handles.Bwareaopen = 1;
handles.Erode_DilateV = 1;
handles.Imfill = 1;
handles.DynamicUpV = 1;

[ handles ] = FeatureExtract( handles.LOGI1,handles.LOGI2, handles );


% hence after crop and imagesegmentation, write it to the file, though this would
% require sharing of pName, fullFileName as handles
dlmwrite(fullfile(handles.pName, handles.predataName), ...
    [handles.cropcor1; handles.cropcor2; handles.ExtractionCor])



%% Verification Section
% Tumourpixelmargin = 25;
% Diahphragmpixelmargin = 20;
% 
% [ handles ] = VerifyExpectation(handles, Tumourpixelmargin, Diahphragmpixelmargin );

%% Plots
handles.SensibleTime = true;

[ handles ] =PlotsCV( handles );


%% Density Plot

    ProbabilityMap = 0;
    for cnt = 1: numel(handles.dicomlist)
        ProbabilityMap = ProbabilityMap + handles.dataMap{cnt};
    end
    % Normalising
    ProbabilityMap = ProbabilityMap/numel(handles.dicomlist);
   
    figure, contourf(ProbabilityMap,'ShowText','on');
axis ij

hTitle  = title ('Probability Density function plot for "patient6 scan1" ');
hXLabel = xlabel('Left-Right (LR) direction '       );
hYLabel = ylabel('Cranial-Caudal (CC) direction'    );


% Adjust font
set( gca                             , 'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel]      ,  'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]               , 'FontSize'   , 10          );
set( hTitle                          , 'FontSize'   , 12          , ...
                                       'FontWeight' , 'bold'      );

 
 %% Scatter Plot, of  handles.TumourYcentroid

figure;
hold on;

% ScatterVMRIY = line((handles.TumourYcentroid)/max(abs(handles.TumourYcentroid)),(handles.DiaphragmYbound)/max(abs(handles.DiaphragmYbound)));
ScatterVMRIY = line((handles.TumourYcentroid/max(abs(handles.TumourYcentroid))),(handles.DiaphragmYbound));
% p1 =     -0.3525  (-0.8803, 0.1752)
% p2 =         6.7  (5.787, 7.614)
% p3 =       2.449  (-0.4117, 5.31)
% x =  -3:0.01:4;
% f(x) = p1*x^2 + p2*x + p3


% Linear model Poly2:
%      f(x) = p1*x^2 + p2*x + p3
% Coefficients (with 95% confidence bounds):
%        p1 =     -0.3525  (-0.8803, 0.1752)
%        p2 =         6.7  (5.787, 7.614)
%        p3 =       2.449  (-0.4117, 5.31)
% 
% Goodness of fit:
%   SSE: 609.1
%   R-square: 0.911
%   Adjusted R-square: 0.9044
%   RMSE: 4.75

% Adjust line properties (functional)
set(ScatterVMRIY                            , ...
  'LineStyle'       , '-.'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );


% Adjust line properties (aesthetics)
set(ScatterVMRIY                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.3 .0 .3]  );


% Add labels
hTitle  = title ('Scatter plot of Diaphragm Belt & MRI Parameter');
hXLabel = xlabel('Normalised CtumourX and CtumourY Parameter\it ( C ) \it '       );
hYLabel = ylabel('Normalised Diaphragm Height Parameter \it ( \mu [normalize])\it'    );


% Add legend
hLegend = legend([ ScatterVMRIY], ...
'Parameters Variation of DiaphragmY and CtumourY in 1 minute'  );

% Adjust font
set( gca                             , 'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel]      ,  'FontName'   , 'AvantGarde');
set([hLegend, gca]                   , 'FontSize'   , 8           );
set([hXLabel, hYLabel]               , 'FontSize'   , 10          );
set( hTitle                          , 'FontSize'   , 12          , ...
                                       'FontWeight' , 'bold'      );

% Adjust axes properties
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'XTick'       , -1:0.5:1.25, ...
  'YTick'       , -1:0.5:1.25, ...
  'LineWidth'   , 1         );
 hold off

%% Scatter Plot, of 

figure;
hold on;
ScatterVMRIX = line((handles.TumourXcentroid)/max(abs(handles.TumourXcentroid)),(handles.DiaphragmYbound)/max(abs(handles.DiaphragmYbound)));

% Adjust line properties (functional)
set(ScatterVMRIX                            , ...
  'LineStyle'       , '-.'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );

% Adjust line properties (aesthetics)
set(ScatterVMRIX                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.0 .7 .7]  );

% Add labels
hTitle  = title ('Scatter plot of Diaphragm Belt & MRI Parameter');
hXLabel = xlabel('Normalised CtumourX and CtumourY Parameter\it ( C ) \it '       );
hYLabel = ylabel('Normalised Diaphragm Height Parameter \it ( \mu [normalize])\it'    );


% Add legend
hLegend = legend([ScatterVMRIX, ], ...
  'Parameters Variation of DiaphragmY and CtumourX in 1 minute' );

% Adjust font
set( gca                             , 'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel]      ,  'FontName'   , 'AvantGarde');
set([hLegend, gca]                   , 'FontSize'   , 8           );
set([hXLabel, hYLabel]               , 'FontSize'   , 10          );
set( hTitle                          , 'FontSize'   , 12          , ...
                                       'FontWeight' , 'bold'      );

% Adjust axes properties
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'XTick'       , -1:0.5:1.25, ...
  'YTick'       , -1:0.5:1.25, ...
  'LineWidth'   , 1         );
 hold off

