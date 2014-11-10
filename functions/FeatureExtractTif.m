function [ tiff ] = FeatureExtract( imageout_t,imageout_d, tiff )
%FEATURE EXTRACTION - this attempt to build the feature extraction as a
%function. The Input parameters are ImgIn that have features to be
%extracted, ExtractionCor, that is the pixel location of the area of
%interest.

% Note, the Cropped Image of Diaphragm should only have one object of interest,

%   This is a modificatino of the feature extraction in GUI.m. Where,
%   handles would pass all relevent arguements, ImgIn is the segmented
%   image to be processed,

plotp = 0;
se = strel('disk',2);
se2 = strel('disk',5);
se3 = strel('disk',20);
bwareaTumour = 30;
bwareaDiaph = 50;


% Check for Exist for substructure, handles.ExtractionCor,
% as normal exist() isn't enough.
doesVarExist = true;% 1 = True
try
    tiff.ExtractionCor;
catch
    doesVarExist = false; % 0 = False
end


for cnt = 1 : tiff.NumberF
    imageout_t{cnt} = imageout_t{cnt}>0; %Cropped Segmentated Image of Tumour
    imageout_d{cnt} = imageout_t{cnt}>0; %Cropped Segmentated Img of Diaphragm
    if cnt == 1
        if doesVarExist == false % If handles.ExtractionsCor doesn't exist
            % Depict the figure for the tumour and diaphragm, to select the
            % ROI of interset to be tracked
            figure, imshow(imageout_t{round(mean([1 tiff.NumberF]))},[])
            [x_t,y_t,vals_t] = impixel;
            close 1 %     %Closing the figure for selection
            figure, imshow(imageout_d{round(mean([1 tiff.NumberF]))},[])
            [x_d,y_d,vals_d] = impixel;
            close 1 %     %Closing the figure for selection
            tiff.ExtractionCor = [x_t y_t x_d y_d];
        else
            % Will be using pre-assigned ExtractionCor
            [x_t,y_t,vals_t] = impixel(imageout_t{cnt},...
                tiff.ExtractionCor(1),tiff.ExtractionCor(2))
            [x_d,y_d,vals_d] = impixel(imageout_d{cnt},...
                tiff.ExtractionCor(3),tiff.ExtractionCor(4))
        end
        
        %generating a blank image size of L
        grain = false(size(imageout_t{cnt}));
        
        A1{cnt}=[x_t,y_t,vals_t(1,:)];
        % Some Image processing operation to fill, remove noise, and
        % artifacts of thresholding
        grain(imageout_t{cnt}==A1{cnt}(3)) = 1;
        % A few if case here to enable the control of checkbox in GUI
        % would be preserve as well I guess
        if tiff.Erode_DilateV == 1;
            grain = imerode(grain, se);
        end
        if tiff.Bwareaopen == 1;
            grain = bwareaopen(grain, bwareaTumour);
        end
        if tiff.Erode_DilateV == 1;
            grain = imdilate(grain, se);
        end
        if tiff.Imfill == 1;
            grain = imfill(grain,'holes');
        end
        dataMap{cnt} = grain;     % passing it to save cropped image with tumour only
        tiff.dataMap{cnt} = grain;
        
        % test for the case if there is more than one non connected threshold region
        % If there is, this would label them differently according to
        % the connetivity, and then the impixel would be able to
        % identify the one selected and isolate it
        grain=bwlabel(grain);
        if max(grain(:))>1 ;
            vals_t=impixel(grain,x_t,y_t);
            dataMap{cnt} = double(grain== vals_t(1));
            tiff.dataMap{cnt} = dataMap{cnt} ;
        end
        
        data=regionprops(tiff.dataMap{cnt},'Area','Centroid','Extrema');
        Atumour(cnt,1)=data(1).Area;
        CtumourX(cnt,1)= data(1).Centroid(1);
        CtumourY(cnt,1)= data(1).Centroid(2);
        Etumour(cnt,1)=data(1).Extrema(2);
        
        %Diaphragm
        imageout_d{cnt} = bwareaopen(imageout_d{cnt},bwareaDiaph);
        imageout_d{cnt} = imdilate(imageout_d{cnt}, se2);
        imageout_d{cnt} = imfill(imageout_d{cnt},'holes');
        imageout_d{cnt} = imerode(imageout_d{cnt}, se2);
        imageout_d{cnt} = imclose(imageout_d{cnt}, se3);
        
        % Finding Boundaries
        [B,L,N] = bwboundaries(imageout_d{cnt});
        coln=find(B{1,1}(:,2)==x_d);
        y_d(cnt) = B{1,1}(coln(1),1);
        ymax  = size(imageout_d{cnt});
                
    else
        %Process for designation of ROI on 2nd run
        % first, establish the location of previous section, and search
        
        if tiff.DynamicUpV == 1;
            try
                [x_t,y_t,vals_t]= impixel(imageout_t{cnt},...
                    CtumourX(cnt-1),(CtumourY(cnt-1)));
            catch exception
                % If it fails to acquire the previous CtumourX, CtumourY, it
                % will try to acquire from the previous previous one.
                fprintf(1, 'Overide at at %d Tumour \n',cnt);
                [x_t,y_t,vals_t]= impixel(imageout_t{cnt},...
                    CtumourX(cnt-2),(CtumourY(cnt-2)));
                continue
            end
        end
        A1{cnt} = [x_t,y_t,vals_t];
    end
    %generating a blank image size of L ( with all 0)
    grain = false(size(imageout_t{cnt}));
    
    
    grain(imageout_t{cnt}==A1{cnt}(3)) = 1;
    if tiff.Erode_DilateV == 1;
        grain = imerode(grain, se);
    end
    if tiff.Bwareaopen == 1;
        grain = bwareaopen(grain, bwareaTumour);
    end
    if tiff.Erode_DilateV == 1;
        grain = imdilate(grain, se);
    end
    if tiff.Imfill == 1;
        grain = imfill(grain,'holes');
    end
    dataMap{cnt} = grain;     % passing it to save cropped image with tumour only
    tiff.dataMap{cnt} = grain;
    % test for the case if there is more than one non connected threshold region
    % If there is, this would label them differently according to
    % the connetivity, and then the impixel would be able to
    % identify the one selected and isolate it
    grain=bwlabel(grain);
    if max(grain(:))>1 ;
        vals_t=impixel(grain,x_t,y_t);
        dataMap{cnt} = double(grain== vals_t(1));
        tiff.dataMap{cnt} = dataMap{cnt} ;
    end
    try              %# Attempt to perform some computation
        data=regionprops(tiff.dataMap{cnt},'Area','Centroid','Extrema');
        Atumour(cnt,1) = data(1).Area;
        CtumourX(cnt,1)= data(1).Centroid(1);
        CtumourY(cnt,1)= data(1).Centroid(2);
        Etumour(cnt,1) = data(1).Extrema(2);
    catch exception  %# Catch the exception
        warning('myfun:warncode','Warning message!')
        fprintf(1, 'Exception at %d Tumour \n',cnt);
        continue       %# Pass control to the next loop iteration
    end
    
    %Diaphragm
    try              %# Attempt to perform some computation
        imageout_d{cnt} = bwareaopen(imageout_d{cnt},bwareaDiaph);
        imageout_d{cnt} = imdilate(imageout_d{cnt}, se2);
        imageout_d{cnt} = imfill(imageout_d{cnt},'holes');
        imageout_d{cnt} = imerode(imageout_d{cnt}, se2);
        imageout_d{cnt} = imclose(imageout_d{cnt}, se3);
        [B,L,N] = bwboundaries(imageout_d{cnt});
        coln=find(B{1,1}(:,2)==x_d);
        y_d(cnt) = B{1,1}(coln(1),1);
        y_max  = size(imageout_d{cnt});
        
    catch exception  %# Catch the exception
        warning('myfun:warncode','Warning message!')
        fprintf(1, 'Exception at %d frame for Diaphragm\n',cnt);
        continue       %# Pass control to the next loop iteration
    end
end



%Changing the extracted data so that it's easier to manage, I should
%perhaps start with these variable from the start, could use some clean
%up in this section.
for cnt = 1 : tiff.NumberF
    cnt
    tiff.Areatumour(cnt) = Atumour(cnt,1);                       %volume of diaphragm
    tiff.CtumourX(cnt) = CtumourX(cnt,1);        %parameters of tumour-x
    tiff.CtumourY(cnt) = CtumourY(cnt,1);       %parameters of tumour-y
    tiff.DiaphragmY(cnt)=y_d(cnt);           %parameters of diaphragm
    tiff.DiaphragmX = x_d;                  %parameters of diaphragm
end

%%%
% hence after crop and imagesegmentation, write it to the file, though this would
% require sharing of pName, fullFileName as handles
dlmwrite(fullfile(tiff.pName, tiff.predataName), ...
    [tiff.cropcor1; tiff.cropcor2; tiff.cropcor3 ; tiff.frame tiff.frame;  tiff.ExtractionCor])



%%% end of features extraction

end

