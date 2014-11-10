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


%% For Images after local Histogram Equalisation

    for cnt = 1 : numel(handles.dicomlist)
        handles.LHI1{cnt} = histeq(handles.I1{cnt});
        handles.LHI2{cnt} = histeq(handles.I2{cnt});
    end
    
 figure, imhist(handles.LHI1{1}), title('after local histeq')
 figure, imshow( handles.LHI1{1},[])

%% For Images after CLAHE

    for cnt = 1 : numel(handles.dicomlist)
        %Apply Contrast-limited Adaptive Histogram Equalization (CLAHE)to
        %an image, ie histeq on [8 8] tiles by default.
        handles.CLdata{cnt} = adapthisteq(handles.data{cnt});
        handles.CLI1{cnt} = imcrop(handles.data{cnt},handles.cropcor1);
        handles.CLI2{cnt} = imcrop(handles.data{cnt},handles.cropcor2);  % times 2 because of the orignal half int
    end
        
figure, imhist(handles.CLI1{1}), title('after CLAHE')
figure, imshow( handles.CLI1{1},[])


% % and the cause for Gaussian filter
%     
%    for cnt = 1 : numel(handles.dicomlist)
%         %Apply Contrast-limited Adaptive Histogram Equalization (CLAHE)to
%         %an image, ie histeq on [8 8] tiles by default.
%         handles.CLdata{cnt} = adapthisteq(handles.data{cnt});
%         handles.CLGI1{cnt} = imcrop(handles.Gdata{cnt},handles.cropcor1);
%         handles.CLGI2{cnt} = imcrop(handles.Gdata{cnt},handles.cropcor2);  % times 2 because of the orignal half int
%    end
   
   %% Gaussian filter
% applying gaussian filter to OriginalLocal histogram and CLAHE processed frame
% just to check the influence 

    for cnt = 1 : numel(handles.dicomlist)
        handles.Gdata{cnt} = Gaussian_fn(handles.data{cnt}, 3,2) ;
        handles.GI1{cnt} = Gaussian_fn(handles.I1{cnt}, 3,2) ;
        handles.GI2{cnt} = Gaussian_fn(handles.I2{cnt}, 3,2) ;
    end  
            
    for cnt = 1 : numel(handles.dicomlist)

        handles.GLHI1{cnt} = Gaussian_fn(handles.LHI1{cnt}, 3,2);
        handles.GLHI2{cnt} = Gaussian_fn(handles.LHI2{cnt}, 3,2);  
    end   

    for cnt = 1 : numel(handles.dicomlist)
        handles.GCLdata{cnt} = Gaussian_fn(handles.CLdata{cnt}, 3,2) ;
        handles.GCLI1{cnt} = Gaussian_fn(handles.CLI1{cnt}, 3,2);
%         handles.GCLI2{cnt} = Gaussian_fn(handles.CLI2{cnt}, 3,2);
    end
    
    
%% Combining Histogram
% As the name suggest, this combine the 1st frame of Original, Local
% Histogram Equalised and CLAHE, compute and have their histogram in one
% graph.  While it might help with human's visual interpretation, 
% I dont think histogram equlisation will have direct influence on
% intensity based segmentation, at least, though the intensity range have
% changed ( seeding, threshold can afford a bigger, lossier grip). However,
% combining this with Gaussian, Gabor filter, (Median filter should be safe)
% might cause some difference, as the range between peaks would be spread 
% to nearby pixels,  thereby affecting the segmentation. 

% The effect/influenced on gradient segmentation would be expected and
% should be examined.

[countsOri,xOri] = imhist(oriI1{1});
% stem(x,counts,'b')
[countsLH,xLH] = imhist(handles.LHI1{1});
% stem(x,counts,'r')
[countsCL,xCL] = imhist(handles.CLI1{1});
% stem(x,counts,'g')
[countsGLH,xGLH] = imhist(handles.GLHI1{1});
[countsGCL,xGCL] = imhist(handles.GCLI1{1});
X = [countsOri, countsLH, countsCL, countsGLH, countsGCL];
Y = [xOri, xLH, xCL, xGLH, xGCL];
figure,stem(Y,X,'MarkerSize',3.75);

% Create legend
legend('Original Img','Local Histogram Img','CLAHE Img', 'LHG', 'CLG') ;
% %Labelling legend
% set(stem1(1),'DisplayName','Original Img');
% set(stem1(2),'DisplayName','Local Histogram Img');
% set(stem1(2),'DisplayName','CLAHE Img');
% Create title
title({'Comparision of Histogram Equalisation Image'});
% Create xlabel
xlabel({'Normalised Bins with respect to range'});
% Create ylabel
ylabel({'Counts'});


%% Test of multiple frame 
% to compare histogram varies between frames
% obviously should be kept within the range of number of frames
% It seems the original histogram stack on top one another too much, hence,
% to compare them, I guess this test would be best to be use on Local
% Histogram, as it has the best spread, however, it is actually only known
% after comparing all of them. Filters such as Gaussian, median is thought 
% to have no influence here.

test = [1 5 13 20 30];

% Hence generate the relevent counts for all test frame
for i = 1 : numel(test);
[y,x] = imhist(handles.LHI1{test(i)});
counts(:,i) = y;
range(:,i) = x;
end

figure,stem(range,counts,'MarkerSize',3.75);
% Create legend
legend('Frame 1','Frame 5','Frame 13','Frame 20','Frame 30') ;
% Create title
title({'Comparision of Different Frame'});
% Create xlabel
xlabel({'Normalised Bins with respect to range'});
% Create ylabel
ylabel({'Counts'});

% Since Frame 5 is obviously the special one here, it would be plot out as
% a figurative compare to others.
for i = 1: 3
figure, imshow( handles.LHI1{test(i)},[])
figure, imshow( handles.data{test(i)},[])
end
% Cause, could potentially be determien by looking at the data ( whole
% dicom image as a whole) as I believe it's cause by the blood circulation
% in heart, despite being surpress....

% Well, apparently, it seems like from the whole dicom image that the
% suppression is caused by blood in aorta arch

%% Computing the number of Outlying peaks on histogram
% Strangely, and perhaps, out of my expectation, on frame 5 there is only 2
% strong peaks, as compare to 3 significant peaks of others. This section
% would check all the frames for number of significant peaks

% Threshold
T = 50;

% Looping through all frames 
for cnt = 1 : numel(handles.dicomlist)
    [y,x] = imhist(handles.LHI1{cnt});
    counts(:,cnt) = y;
    range(:,cnt) = x;
    Map =  counts > T;
    PeakCnt(cnt) = sum(Map(:,cnt));
end

%Plotting
createbarchart(PeakCnt)
% % Create title
% title({'Number of peaks in different frame'});
% % Create xlabel
% xlabel({'Frame'});
% % Create ylabel
% ylabel({'Counts of significant peaks'});

% Therefore, we can conclude that histogram normalisation, or perhaps, more
% accurately, histogram matching is required. Though, due to the chances of
% phase, matching all image to the very first frame isn't very adequate

% %%% Let the figures dance!!!
% set(0,'DefaultFigureVisible','on')

%% Segmentation of Tumour

%% Segmentation method 1: Global Otsu method
% testing waters, though there are a little bit too numerous methods
% by my understanding, due to the distribution different oh peaks on globally
% and locally, there would be a difference.

handles.OtsuBinV = 3;

for cnt = 1 : numel(handles.dicomlist)
    % Acting on Original data
    handles.GOdataT{cnt} = otsu(handles.data{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    handles.GOI1{cnt} = imcrop(handles.GOdataT{cnt},handles.cropcor1)*handles.OtsuBinV;
    % Acting on Gaussian filtered data
    handles.GOGdataT{cnt} = otsu(handles.Gdata{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    handles.GOGI1{cnt} = imcrop(handles.GOGdataT{cnt},handles.cropcor1)*handles.OtsuBinV;
    % Acting on CLAHE data
    handles.GOCLdataT{cnt} = otsu(handles.CLdata{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    handles.GOCLI1{cnt} = imcrop(handles.GOCLdataT{cnt},handles.cropcor1)*handles.OtsuBinV;
    % Acting on  Gaussian CLAHE data
    handles.GOGCLdataT{cnt} = otsu(handles.GCLdata{cnt},handles.OtsuBinV)*handles.OtsuBinV;
    handles.GOGCLI1{cnt} = imcrop(handles.GOGCLdataT{cnt},handles.cropcor1)*handles.OtsuBinV;
    %     handles.I2{cnt} = imcrop(handles.dataT{cnt},handles.cropcor2)*handles.OtsuBinV;  % times 2 because of the orignal half int
end

% going through to depict imshow on test frames, 
% too many figures is bad bad bad.
for i = 1 : 3
    figure('Name',' Original'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
    figure('Name','Global Otsu method on Original'),...
        imshow(handles.GOI1{test(i)},[],'InitialMagnification',100), title({'Global Otsu method on Original'})
    figure('Name','Global Otsu method on Gaussian data'),...
        imshow(handles.GOGI1{test(i)},[],'InitialMagnification',100), title({'Global Otsu method on Gaussian data'})
    figure('Name','Global Otsu method on CLAHE data'),...
        imshow(handles.GOCLI1{test(i)},[],'InitialMagnification',100), title({'Global Otsu method on CLAHE data'})
    figure('Name','Global Otsu method on Gaussian CLAHE data'),...
        imshow(handles.GOGCLI1{test(i)},[],'InitialMagnification',100), title({'Global Otsu method on Gaussian CLAHE data'})
end

% This depicted the importance of control on Otsu bin, the accuracy and
% influence of Gaussian, CLAHE method etc. after CLAHE segmentation is much
% better

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
for i = 1 : 3
%     figure('Name',' Original'),...
%         imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
%     figure('Name','Local Otsu method on Original'),...
%         imshow(handles.LOI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Original'}),
    figure('Name','Local Otsu method on Gaussian data'),...
        imshow(handles.LOGI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian data'})
%     figure('Name','Local Otsu method on LH data'),...
%         imshow(handles.LOLHI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on LH data'})
%     figure('Name','Local Otsu method on Gaussian LH data'),...
%         imshow(handles.LOGLHI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian LH data'})
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

for i = 1: 3
    figure('Name',' Original'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
    figure('Name','MCWS on Gaussian data'),...
        imshow(handles.MCWSI1{test(i)},[],'InitialMagnification',100), title({'MCWS on Gaussian data'})
%        figure('Name','MCWS on   data'),...
%         imshow(handles.MCWSI1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian LH data'})
end

% Gaussian filter would bode bad for Gradient chances, would perhaps affect
% MCWS, too troublesome to examine now, would have to make region merging into
% another function

%% Segmentation method 4:Template Matching method


% Loading template
template = double(imread([handles.ImageFolder 'template.png']));
template = template(:,:,1);
%          template=rgb2gray(double(template>0));

for cnt = 1 : numel(handles.dicomlist)
    [BestRow, BestCol] = MaxCorrelation2(handles.I1{cnt}, template);
    % inserting template at the correct location of the coronal plane
    % of picture
    [hT,wT] = size(template);
    handles.TMI1{cnt} = false(size(handles.I1{1}));
    handles.TMI1{cnt}(BestRow : BestRow + hT - 1, BestCol : BestCol + wT- 1 ) = template > 0;
end

%Plot of Original Frame with template overlay on top of it's frame
for i = 1: 3
    figure('Name',' Original'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({' Original'})
    
    green = cat(3, zeros(size(handles.I1{1})), ones(size(handles.I1{1})), zeros(size(handles.I1{1})));
    hold on
    h = imshow(green);
    % Use  TM1 map as the AlphaData for the solid green image overlay on
    % Original data
    set(h, 'AlphaData', handles.TMI1{test(i)}*0.3)
    hold off
    %     figure('Name','Template Matching method on Original'),...
    %         imshow(handles.TMI1{test(i)},[],'InitialMagnification',100), title({'Template Matching method on Original'})
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

[ handles ] = FeatureExtract( handles.LOGI1,handles.LOI2, handles );


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
