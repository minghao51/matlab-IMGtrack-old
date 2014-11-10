function [ handles ] = FeatureExtract( ImgIn_t,ImgIn_d, handles )
%FEATURE EXTRACTION - this attempt to build the feature extraction as a
%function. The Input parameters are ImgIn that have features to be
%extracted, ExtractionCor, that is the pixel location of the area of
%interest.

% % % % % old Parameters, moving these sections to into the loops
% % % % % for cnt = 1: numel(handles.dicomlist)
% % % % %     handles.Areatumour(cnt) = Atumour(cnt,1);                        %volume of diaphragm
% % % % %     handles.CtumourX(cnt) = CtumourX(cnt,1);        %parameters of tumour-x
% % % % %     handles.CtumourY(cnt) = CtumourY(cnt,i);        %parameters of tumour-y
% % % % %     handles.DiaphragmY(cnt)=y_d(cnt);           %parameters of diaphragm
% % % % %     handles.DiaphragmX = x_d;                  %parameters of diaphragm
% % % % % end


% Note, the Cropped Image of Diaphragm should only have one object of interest,

%   This is a modificatino of the feature extraction in GUI.m. Where,
%   handles would pass all relevent arguements, ImgIn is the segmented
%   image to be processed,

N_ROI = 1;
plotp = 0;
se = strel('disk',2);
se2 = strel('disk',15);
% se2/3 = strel('disk',5); for patient6scan1, 15 for corpat case
se3 = se2;
bwareaTumour = 30;
bwareaDiaph = 10;
%600 for tumour, 200 for RESP case, 10 for corpat
% 30 for shaba case
% 150 for all others causes

% Optimising
%prelocating dataMap, the final extracted map of tumour
dataMap{numel(handles.dicomlist)} =  {};
handles.dataMap{numel(handles.dicomlist)} =  {};
handles.dataMapDiaphragm{numel(handles.dicomlist)} =  {};

% Check for Exist for substructure, handles.ExtractionCor,
% as normal exist() isn't enough.
doesVarExist = true;% 1 = True
try
    handles.ExtractionCor;
catch
    doesVarExist = false; % 0 = False
end


for cnt = 1 : numel(handles.dicomlist)
    fprintf(1, 'intitialing at %d Tumour \n',cnt);
    imageout_t{cnt} = ImgIn_t{cnt}; %Cropped Segmentated Image of Tumour
    imageout_d{cnt} = ImgIn_d{cnt}; %Cropped Segmentated Img of Diaphragm
    if cnt == 1
        if doesVarExist == false % If handles.ExtractionsCor doesn't exist
            % Depict the figure for the tumour and diaphragm, to select the
            % ROI of interset to be tracked
            figure, imshow(imageout_t{cnt},[])
            [x_t,y_t,vals_t] = impixel;
            close 1 %     %Closing the figure for selection
            figure, imshow(imageout_d{cnt},[])
            [x_d,y_d,vals_d] = impixel;
            close 1 %     %Closing the figure for selection
            handles.ExtractionCor = [x_t y_t x_d y_d];
        else
            % Will be using pre-assigned ExtractionCor
            [x_t,y_t,vals_t] = impixel(imageout_t{cnt},...
                handles.ExtractionCor(1),handles.ExtractionCor(2));
            [x_d,y_d,vals_d] = impixel(imageout_d{cnt},...
                handles.ExtractionCor(3),handles.ExtractionCor(4));
        end
        
        
        %loop to show ROI only and calculate their area
        %generating a blank image size of L
        grain = false(size(imageout_t{cnt}));
        
        for i = 1:N_ROI
            A1{cnt,i}=[x_t(i),y_t(i),vals_t(i,:)];
            % Some Image processing operation to fill, remove noise, and
            % artifacts of thresholding
            grain(imageout_t{cnt}==A1{cnt,i}(3)) = 1;
            % A few if case here to enable the control of checkbox in GUI
            % would be preserve as well I guess
            if handles.Erode_DilateV == 1;
                grain = imerode(grain, se);
            end
            if handles.Bwareaopen == 1;
                grain = bwareaopen(grain, bwareaTumour);
            end
            if handles.Erode_DilateV == 1;
                grain = imdilate(grain, se);
            end
            if handles.Imfill == 1;
                grain = imfill(grain,'holes');
            end
            dataMap{cnt} = grain;     % passing it to save cropped image with tumour only
            handles.dataMap{cnt} = grain;
            
            % test for the case if there is more than one non connected threshold region
            % If there is, this would label them differently according to
            % the connetivity, and then the impixel would be able to
            % identify the one selected and isolate it
            grain=bwlabel(grain);
            if max(grain(:))>1 ;
                vals_t=impixel(grain,x_t,y_t);
                dataMap{cnt} = double(grain== vals_t(1));
                handles.dataMap{cnt} = dataMap{cnt} ;
            end
            
            data=regionprops(handles.dataMap{cnt},'Area','Centroid','Extrema');
            fprintf(1, 'regionprops at %d Tumour \n',cnt);
            Atumour(cnt,i)=data(1).Area;
            CtumourX(cnt,i)= data(1).Centroid(1);
            CtumourY(cnt,i)= data(1).Centroid(2);
            Etumour(cnt,i)=data(1).Extrema(2);
            
            % passing the parameters to be shared between routines
            handles.Areatumour(cnt) = Atumour(cnt,i);
            handles.CtumourX(cnt)= CtumourX(cnt,i); 
            handles.CtumourY(cnt)= CtumourY(cnt,i); 
            
            figure, imshow(grain), title(' Centroid Locations')
            hold on
            for l = 1 : N_ROI
                plot(x_t(l),y_t(l),'ro') % Plotting impixel
                plot(data(l).Centroid(1), data(l).Centroid(2),'bo');
            end
            hold off
            
            
            %             % In case that the centroid doens't clip to the Region of
            %             % Interest, Here the CtumourX, CtumourY, Ediaphragm of all is
            %             % set to the 1st coordinates, this only works if the ROI
            %             % doesn't move much from its initial coordiantes
            %             for cnt = 1: numel(handles.dicomlist)
            %                 CtumourX(cnt,i)= data(vals(i,1)).Centroid(1);
            %                 CtumourY(cnt,i)= data(vals(i,1)).Centroid(2);
            %                 Ediaphragm(cnt,i)=data(vals(i,1)).Extrema(2);
            %             end
        end
        
        %
        %         % Plotting the Grain( ImageSegmentation-selected) ROI with centroid marked on
        %         % figure for Tumour
        %         if plotp == 2
        %             figure, imshow(grain), title(' Centroid Locations')
        %             hold on
        %             for l = 1 : N_ROI
        %                 plot(data(l).Centroid(1), data(l).Centroid(2),'bo');
        %             end
        %             hold off
        %         end
        
        % Plotting the Grain( ImageSegmentation-selected) ROI with centroid marked on
        %         % figure for Diaphragm
        %         if plotp == 2
        %             figure,imshow(imageout_d{cnt},[])
        %         end
        % Removing noise
        imageout_d{cnt} = bwareaopen(imageout_d{cnt},bwareaDiaph);
        imageout_d{cnt} = imdilate(imageout_d{cnt}, se2);
        imageout_d{cnt} = imfill(imageout_d{cnt},'holes');
        imageout_d{cnt} = imerode(imageout_d{cnt}, se2);
        imageout_d{cnt} = imclose(imageout_d{cnt}, se3);
        %         imageout_d{cnt} = imclose(imageout_d{cnt}, se3);
        handles.dataMapDiaphragm{cnt}= imageout_d{cnt};
        %         imageout_d{cnt} = bwlabel(imageout_d{cnt}, se3);
        % Finding Boundaries
        [B,L,N] = bwboundaries(imageout_d{cnt});
        %changed B{1,1} to B{2,1}
        coln=find(B{1,1}(:,2)==x_d);
        y_d(cnt) = B{1,1}(coln(1),1);
        ymax  = size(imageout_d{cnt});
        
        %         hold on
        %         line([x_d,x_d],[0,ymax(1)],'Color','r','LineWidth',2)
        %         plot(x_d,y_d(cnt),'bo');
        %         hold off
        
        
    else
        %Process for designation of ROI on 2nd run
        % first, establish the location of previous section, and search
        for i = 1:N_ROI
            if handles.DynamicUpV == 1;
                try
                    [x_t,y_t,vals_t]= impixel(imageout_t{cnt},...
                        CtumourX(cnt-1,i),(CtumourY(cnt-1,i)));
                    %                     fprintf(1, 'Assignemnt of Impixel at %d Tumour \n',cnt);
                catch exception
                    % If it fails to acquire the previous CtumourX, CtumourY, it
                    % will try to acquire from the previous previous one.
                    fprintf(1, 'Overide at at %d Tumour \n',cnt);
                    [x_t,y_t,vals_t]= impixel(imageout_t{cnt},...
                    CtumourX(cnt-3,i),(CtumourY(cnt-3,i)));
                %                         CtumourX(cnt-2,i),(CtumourY(cnt-2,i)));
                end
            end
            A1{cnt,i} = [x_t,y_t,vals_t];
        end
        %generating a blank image size of L ( with all 0)
        grain = false(size(imageout_t{cnt}));
        
        %loop to show ROI only and calculate their area
        for i = 1:N_ROI
            %set the ROI interest to be 1
            grain(imageout_t{cnt}==A1{cnt,i}(3)) = 1;
            if handles.Erode_DilateV == 1;
                grain = imerode(grain, se);
            end
            if handles.Bwareaopen == 1;
                grain = bwareaopen(grain, bwareaTumour);
            end
            if handles.Erode_DilateV == 1;
                grain = imdilate(grain, se);
            end
            if handles.Imfill == 1;
                grain = imfill(grain,'holes');
            end
            dataMap{cnt} = grain;     % passing it to save cropped image with tumour only
            handles.dataMap{cnt} = grain;
            % test for the case if there is more than one non connected threshold region
            % If there is, this would label them differently according to
            % the connetivity, and then the impixel would be able to
            % identify the one selected and isolate it
            grain=bwlabel(grain);
            if max(grain(:))>1 ;
                vals_t=impixel(grain,x_t,y_t);
                dataMap{cnt} = double(grain== vals_t(1));
                handles.dataMap{cnt} = dataMap{cnt} ;
            end
            try              %# Attempt to perform some computation
                data=regionprops(handles.dataMap{cnt},'Area','Centroid','Extrema');
                Atumour(cnt,i) = data(1).Area;
                CtumourX(cnt,i)= data(1).Centroid(1);
                CtumourY(cnt,i)= data(1).Centroid(2);
                Etumour(cnt,i) = data(1).Extrema(2);
                %                 fprintf(1, 'Pass at %d Tumour \n',cnt);
                
                
                % passing the parameters to be shared between routines
                handles.Areatumour(cnt) = Atumour(cnt,i) ;
                handles.CtumourX(cnt)= CtumourX(cnt,i);
                handles.CtumourY(cnt)= CtumourY(cnt,i);
                
            catch exception  %# Catch the exception
                warning('myfun:warncode','Warning message!')
                fprintf(1, 'Exception at %d Tumour \n',cnt);
                continue       %# Pass control to the next loop iteration
            end
        end
        
        % Plot the Grain( ImageSegmentation-selected) ROI with centroid marked on
        % figure for Tumour
        %             if plotp == 2
        %                  figure, imshow(grain), title(' Centroid Locations')
        %                 hold on
        %                 for l = 1 : N_ROI
        %                     plot(x_t(l),y_t(l),'ro') % Plotting impixel
        %                     plot(data(l).Centroid(1), data(l).Centroid(2),'bo');
        %                 end
        %                 hold off
        %             end
        
        %             % Plotting the Grain( ImageSegmentation-selected) ROI with centroid marked on
        %             % figure for Diaphragm
        %             if plotp == 2
        %                 figure,imshow(imageout_d{cnt},[])
        %             end
        
        try              %# Attempt to perform some computation
            imageout_d{cnt} = bwareaopen(imageout_d{cnt},bwareaDiaph);
            imageout_d{cnt} = imdilate(imageout_d{cnt}, se2);
            imageout_d{cnt} = imfill(imageout_d{cnt},'holes');
            imageout_d{cnt} = imerode(imageout_d{cnt}, se2);
            imageout_d{cnt} = imclose(imageout_d{cnt}, se3);
            handles.dataMapDiaphragm{cnt}= imageout_d{cnt};
            [B,L,N] = bwboundaries(imageout_d{cnt});
            coln=find(B{1,1}(:,2)==x_d);
            y_d(cnt) = B{1,1}(coln(1),1);
            y_max  = size(imageout_d{cnt});
            %                 This would plot the line, and the highest point of
            %                 diaphgragm on every fram
            %                 if plotp == 2
            %                     hold on
            %                     line([x_d,x_d],[0,ymax(1)],'Color','r','LineWidth',2)
            %                     plot(x_d,y_d(cnt),'bo');
            %                     hold off
            %                 end
            
            %
            handles.DiaphragmY(cnt)=y_d(cnt);
            handles.DiaphragmX = x_d;
            
        catch exception  %# Catch the exception
            warning('myfun:warncode','Warning message!')
            fprintf(1, 'Exception at %d frame for Diaphragm\n',cnt);
            continue       %# Pass control to the next loop iteration
        end
    end
end


%Changing the extracted data so that it's easier to manage, I should
%perhaps start with these variable from the start, could use some clean
%up in this section.

% % % % % old Parameters, moving these sections to above
% % % % % for cnt = 1: numel(handles.dicomlist)
% % % % %     handles.Areatumour(cnt) = Atumour(cnt,1);                        %volume of diaphragm
% % % % %     handles.CtumourX(cnt) = CtumourX(cnt,1);        %parameters of tumour-x
% % % % %     handles.CtumourY(cnt) = CtumourY(cnt,i);        %parameters of tumour-y
% % % % %     handles.DiaphragmY(cnt)=y_d(cnt);           %parameters of diaphragm
% % % % %     handles.DiaphragmX = x_d;                  %parameters of diaphragm
% % % % % end


% debuging check
fprintf(1, ' %d Tumour \n',numel(handles.CtumourX));

fprintf(1, ' %d Diaphragm \n',numel(handles.DiaphragmY));

%%%
% hence after crop and imagesegmentation, write it to the file, though this would
% require sharing of pName, fullFileName as handles
dlmwrite(fullfile(handles.pName, handles.predataName), ...
    [handles.cropcor1; handles.cropcor2; handles.ExtractionCor])



%%% end of features extraction

end

