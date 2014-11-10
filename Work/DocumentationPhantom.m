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

% loop to repeat the loading of different phantom sets
for phantomsets=1:4
    
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
            
%% Test of multiple frame 
test = [1 5 13 20 30];
    
   

%% Segmentation of Tumour
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
    %         handles.I1{cnt} = imcrop(handles.dataT{cnt},handles.cropcor1)*2;
    %         handles.I2{cnt} = imcrop(handles.dataT{cnt},handles.cropcor2)*2;  % times 2 because of the orignal half int
    % Acting on Original data
    handles.LOI1{cnt} = otsu(handles.I1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    % Acting on Gaussian filtered data
    handles.LOGI1{cnt} = otsu(handles.GI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    % Acting on LH data
%     handles.LOLHI1{cnt} = otsu(handles.LHI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    % Acting on LHG data
%     handles.LOGLHI1{cnt} = otsu(handles.GLHI1{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    %     handles.I2{cnt} = otsu(handles.I2{cnt},handles.OtsuBinV)*handles.OtsuBinV;
end

% going through to depict imshow on test frames, 
% for i = 1 : 3
%     figure('Name',' Original'),...
%         imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
%     figure('Name','Local Otsu method on Original'),...
%         imshow(handles.LOI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Original'}),
%     figure('Name','Local Otsu method on Gaussian data'),...
%         imshow(handles.LOGI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian data'})
% end

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
% Since currently plotting by acquisition time and content time on the siemens MRI
% makes little sense, the time plot is just, far from resonable
handles.SensibleTime = false;
% The number of frames starting from the first that which makes sense in
% its temporal time.
handles.TimeSensibleRange = 2;

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
end
