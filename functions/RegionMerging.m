function [ handles ] = RegionMerging( handles )
%REGION MERGING - this will use regionprops to measure the Mean and
%Variance of each segmentated region, and determine the Mean and Max
%intensity, allowing region within varience of the greater regino to merge
%to farm a uniform region, this is done to counteract the oversegmentation
%of MCWS method.
%   Detailed explanation goes here


%Region merging, loop over all frame (cnt) and then small regions(k).
for cnt = 1 : numel(handles.dicomlist)
    
    % Acquire the pixelvalues etc using region props
    RegionMeasurement = regionprops(handles.MCWSI1{cnt},  handles.I1{cnt}, 'PixelValues',...
        'MeanIntensity', 'MaxIntensity', 'MinIntensity', 'PixelIdxList' );
    numberOfRegions(cnt) = size(RegionMeasurement, 1);
    
    % Asign the extracted values, and calculate the mean, std etc
    for k = 1 : numberOfRegions(cnt)
        Pixelmean( k,cnt) = RegionMeasurement(k).MeanIntensity;
        PixelMaxValues( k,cnt) = RegionMeasurement(k).MaxIntensity;
        PixelMinValues( k,cnt) = RegionMeasurement(k).MinIntensity;
        PixelVar( k,cnt) = var(double(RegionMeasurement(k).PixelValues));
        Pixelstd( k,cnt) = std(double(RegionMeasurement(k).PixelValues));
    end
end


% Loop over all frame and region, determine if the region is worth
% merging based on the mean and var value.
for cnt = 1 : numel(handles.dicomlist) %
    % Loop over all region, comparing to themself one by one
    exemption (1:numberOfRegions(cnt)) = 0;
    mergedMask= 0;
    for k = 1 : numberOfRegions(cnt)
        if k ~= exemption
            mergedMask = handles.MCWSI1{cnt} == k;
            for j = 2 : numberOfRegions(cnt)
                if j ~= exemption
                    % if they are to be merged
                    if abs(Pixelmean(j, cnt) - Pixelmean(k, cnt)) < max(Pixelstd(k, cnt)+ Pixelstd(j, cnt))
                        %                         [abs(Pixelmean(j, cnt) - Pixelmean(k, cnt)) Pixelstd(k, cnt)]
                        %                    mergedMask = region1Mask & region2Mask; % Assumes binary images.
                        mergedMask = int8(mergedMask) + int8(handles.MCWSI1{cnt} == j); % Assumes binary images.
                        %                             figure, imshow (mergedMask)
                        % the j value should be exempted from k
                        exemption (j) = j;
                    end
                    handles.MCWSI1{cnt}(mergedMask>0) = k;
                end
                %                     mergedImage = handles.MCWSI1{cnt}; % Initialize new image.
                %                     mergedImage(~mergedMask) = 0; % Zero out everything outside the merged region.
                %                     handles.MCWSI1{cnt}(mergedMask) = k;
                %                 figure, imshow(handles.MCWSI1{cnt},[])
            end
        end
    end
    %     figure, imshow(handles.MCWSI1{cnt},[]), title('before median filter')
%     handles.MCWSI1{cnt} = medfilt2(handles.MCWSI1{cnt},[5 5]);
%     figure, imshow(handles.MCWSI1{cnt},[]), title('after median filter')
    %     homerunc(cnt) = bwlabel(handles.MCWSI1{cnt},8);
    %     figure, imshow(homerunc(cnt),[])
end

% bwlabel(handles.MCWSI1{cnt},8)

end

