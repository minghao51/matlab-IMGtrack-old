%Otsu
for i = 1 : 7
    h=figure('Name','Local Otsu method on Gaussian data'),...
        imshow(handles.I1{test(i)},[],'InitialMagnification',100), title({'Local Otsu method on Gaussian data'});
    [BOtsu,LOtsu] = bwboundaries(handles.LOGI1{test(i)},'noholes');
    [BMaxC,LMaxC] = bwboundaries(handles.TMI1{test(i)},'noholes');
    [BSeed,LSeed] = bwboundaries(handles.RSI1{test(i)},'noholes');
    [BMCWS,LMCWS] = bwboundaries(handles.MCWSI1{test(i)},'noholes');
    hold on
%     for k = 1:length(BOtsu)
%         boundaryOtsu = BOtsu{k};
%         Otsu = plot(boundaryOtsu(:,2), boundaryOtsu(:,1), 'r', 'LineWidth', LineWidth);
%     end
%     for k = 1:length(BMaxC)
%         boundaryMaxC = BMaxC{k};
%         MaxC= plot(boundaryMaxC(:,2), boundaryMaxC(:,1), 'y', 'LineWidth', LineWidth);
%     end
%     for k = 1:length(BSeed)
%         boundarySeed = BSeed{k};
%         Seed = plot(boundarySeed(:,2), boundarySeed(:,1), 'b', 'LineWidth', LineWidth);
%     end
     for k = 1:length(BMCWS)
        boundaryMCWS = BMCWS{k};
        MCWS = plot(boundaryMCWS(:,2), boundaryMCWS(:,1), 'g', 'LineWidth', LineWidth);
    end
    hLegend = legend([MCWS ], ...%Otsu,MaxC,Seed,
'Watershed Segmentation' );%' Otsu segmentation', 'TemplateSegmentation' ,'Seed Segmentation',

    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
   
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\CombinedSeg\mcws\' handles.dicomlist(i).name '.png']) ;
end

%MCWS
for i = 1: 7
    L =MCWS1(handles.GI1{test(i)}, bgmV, MBGM, fgmV, MFGM, imregionalmaxV );
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    h = figure, imshow(Lrgb,[]);
    title('Colored watershed label matrix (Lrgb)')
    
    saveas(h,['C:\Users\Minghao\Desktop\research project\Thesis\MCWSSeg\' handles.dicomlist(i).name '.png'])
end

%MaxCorr
  
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

% Seedgrowing
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