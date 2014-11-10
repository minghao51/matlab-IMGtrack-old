%% Documents
%Test case for programs in GUI: loading all files as usual, but stop
% before segmenting, would be changes to test all different segmentation
% methods

%% Clearing and some default settings
clear all
close all

% docked figures are easier to manage
% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')

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

% %Depiciting the Image
% figure, imhist(handles.I1{1}), title('Original Image Histogram')
% figure, imshow( handles.I1{1},[])
   
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
    handles.LOGI1Bin2{cnt} = otsu(handles.GI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
end

% handles.OtsuBinV = 3;
% for cnt = 1 : numel(handles.dicomlist)
%     handles.LOGI1Bin3{cnt} = otsu(handles.GI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
% end

figure,imshow(handles.GI1{1},[])
figure,imshow(handles.LOGI1Bin2{1},[])
% figure,imshow(handles.LOGI1Bin3{2},[])


% % going through to depict imshow on test frames, 
% for i = 1 : 7
%     h=figure('Name','Local Otsu method on Gaussian data'),...
%         imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian data'});
%     [B,L] = bwboundaries(handles.LOGI1{test(i)},'noholes');
%     hold on
%     for k = 1:length(B)
%         boundary = B{k};
%         plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', LineWidth)
%     end
%     hold off
%     set(gca,'position',[0 0 1 1],'units','normalized')
%    
%     saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\OtsuSegmentation\' handles.dicomlist(i).name '.png']) ;
% end



%% Segmentation of Diaphragm

handles.OtsuBinV = 2;

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

[ handles ] = FeatureExtract( handles.LOGI1Bin2,handles.LOGI2, handles );


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
handles.SensibleTime = false;
handles.TimeSensibleRange = 2;

    [ handles ] =PlotsCV( handles );

%%

[r,p]=corrcoef(handles.TumourXcentroid,handles.DiaphragmYbound)
CtumourX=handles.TumourXcentroid;
CtumourY=handles.TumourYcentroid;
DiaphragmY=handles.DiaphragmYbound;
Time= handles.time_axis;

%% Density Plot PDF

    ProbabilityMap = 0;
    for cnt = 1: numel(handles.dicomlist)
        ProbabilityMap = ProbabilityMap + handles.dataMap{cnt};
    end
    % Normalising
    ProbabilityMap = ProbabilityMap/numel(handles.dicomlist);
   
    figure, contourf(ProbabilityMap,'ShowText','on');
axis ij

hTitle  = title ('Probability Density function plot for "Khamis-1-Sag" ');
hXLabel = xlabel('(AP) direction '       );
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
ScatterVMRIY = line((handles.TumourYcentroid),(handles.DiaphragmYbound));
% ScatterVMRIY_period1 = line(handles.TumourYcentroid(2:5),(handles.DiaphragmYbound(2:5)));
% ScatterVMRIY_period2 = line(handles.TumourYcentroid(5:9),(handles.DiaphragmYbound(5:9)));
% ScatterVMRIY_period3 = line(handles.TumourYcentroid(9:11),(handles.DiaphragmYbound(9:11)));

% Adjust line properties (functional)
set(ScatterVMRIY                            , ...
  'LineStyle'       , '-.'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );
% set(ScatterVMRIY_period1                          , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [0 0 .5]    );
% set(ScatterVMRIY_period2                           , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [0 .5 0]    );
% set(ScatterVMRIY_period3                           , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [.5 0 0]    );

% Adjust line properties (aesthetics)
set(ScatterVMRIY                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.3 .0 .3]  );
% set(ScatterVMRIY_period1                          , ...
%   'LineWidth'       , 2           );
% set(ScatterVMRIY_period2                          , ...
%   'LineWidth'       , 2           );
% set(ScatterVMRIY_period3                          , ...
%   'LineWidth'       , 2           );

% Add labels
hTitle  = title ('Scatter plot of CtumourY & Diaphragm Height');
hXLabel = xlabel(' CtumourY '       );
hYLabel = ylabel(' Diaphragm Height '   );

% Add legend
hLegend = legend([ScatterVMRIY], ...
  'Parameters Varixation of DiaphragmY and CtumourY in 1 minute'  );

% hLegend = legend([ScatterVMRIX, ScatterVMRIY_period1, ScatterVMRIY_period2, ScatterVMRIY_period3 ], ...
%   'Parameters Varixation of DiaphragmY and CtumourY in 1 minute' , ...
%   ' 1st respiration period' , ...
%   '2nd respiration period', ...
%   '3rd respiration period' );

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
  'LineWidth'   , 1         );
 hold off

%% Scatter Plot, of 

figure;
hold on;
ScatterVMRIX = line((handles.TumourXcentroid),(handles.DiaphragmYbound));
% ScatterVMRIY_period1 = line(handles.TumourXcentroid(2:5),(handles.DiaphragmYbound(2:5)));
% ScatterVMRIY_period2 = line(handles.TumourXcentroid(5:9),(handles.DiaphragmYbound(5:9)));
% ScatterVMRIY_period3 = line(handles.TumourXcentroid(9:14),(handles.DiaphragmYbound(9:14)));

% Adjust line properties (functional)
set(ScatterVMRIX                            , ...
  'LineStyle'       , '-.'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );
% set(ScatterVMRIY_period1                          , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [0 0 .5]    );
% set(ScatterVMRIY_period2                           , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [0 .5 0]    );
% set(ScatterVMRIY_period3                           , ...
%   'Marker'          , '+'         , ...
%   'Color'           , [.5 0 0]    );

% Adjust line properties (aesthetics)
set(ScatterVMRIX                            , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.0 .7 .7]  );
% set(ScatterVMRIY_period1                          , ...
%   'LineWidth'       , 2           );
% set(ScatterVMRIY_period2                          , ...
%   'LineWidth'       , 2           );
% set(ScatterVMRIY_period3                          , ...
%   'LineWidth'       , 2           );

% Add labels
hTitle  = title ('Scatter plot of CtumourX & Diaphragm Height');
hXLabel = xlabel(' CtumourX  '       );
hYLabel = ylabel(' Diaphragm Height '    );

% Add legend
hLegend = legend([ScatterVMRIX ], ...
  'Parameters Varixation of DiaphragmY and CtumourX in 1 minute' );

% hLegend = legend([ScatterVMRIX, ScatterVMRIY_period1, ScatterVMRIY_period2, ScatterVMRIY_period3 ], ...
%   'Parameters Varixation of DiaphragmY and CtumourX in 1 minute' , ...
%   ' 1st respiration period' ,...
%   '2nd respiration period',...
%   '3rd respiration period' );

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
  'LineWidth'   , 1         );
 hold off

