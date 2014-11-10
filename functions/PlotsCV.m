function [ handles ] = PlotsCV( handles )
%% Plotting Centroid position changes, diphragm changes and Volume
% This section establish the actual displacement and time parameters as
% well, instead of just relating to pixel coordinates and time.

% if(exist('handles.SensibleTime','var')==0), handles.SensibleTime=true; end

% %% Prelocating x_axis, time
% x_axis{numel(handles.dicomlist)} =  {};


if handles.SensibleTime == true;
    % info.AcquisitionTime
    for cnt = 1 : numel(handles.dicomlist)
        %     InitialTime =  handles.info{1}.AcquisitionTime;
        x_axis{cnt} = str2num(handles.info{cnt}.AcquisitionTime ) - str2num(handles.info{1}.AcquisitionTime);
        %  x1{cnt} = x1{cnt}
    end
    handles.time_axis  = cell2mat (x_axis);
    
else
    
    for cnt = 1 : numel(handles.dicomlist)
        %     InitialTime =  handles.info{1}.AcquisitionTime;
        x_axis{cnt} = (str2num(handles.info{cnt}.AcquisitionTime ) - str2num(handles.info{1}.AcquisitionTime));
        %  x1{cnt} = x1{cnt}
    end
    temporaltime = (x_axis{handles.TimeSensibleRange}-x_axis{1})/(handles.TimeSensibleRange-1);
    endpoint = (numel(handles.dicomlist)-1)*temporaltime;
%     handles.time_axis  = linspace(0,endpoint-1,numel(handles.dicomlist));
% handles.time_axis  = linspace(0,endpoint,numel(handles.dicomlist)+1);
handles.time_axis = (0:temporaltime:endpoint);
end

%% ACtual physical Displacement
%transforming the centroid pixels locations to actual displacement and time

%info.pixelspacing
handles.TumourXcentroid = handles.CtumourX;
handles.TumourXcentroid = (handles.TumourXcentroid - nanmean(handles.CtumourX))*handles.info{1}.PixelSpacing(1); %Mean ignoring NaN values
handles.TumourYcentroid = handles.CtumourY;
handles.TumourYcentroid = (handles.TumourYcentroid - nanmean(handles.CtumourY))*handles.info{1}.PixelSpacing(1);
handles.DiaphragmYbound = handles.DiaphragmY;
handles.DiaphragmYbound = (handles.DiaphragmYbound - nanmean(handles.DiaphragmY))*handles.info{1}.PixelSpacing(1);

%% Plot of Centroid of tumour and Diaphragm coordinates
% This plot the appropriate data, with the respectice time in s and
% displacement in mm

figure
hold on
% hleg1 = legend('Ctumour-x AP/RL','Ctumour-y CC','Diaphragm-y CC');
TumourXcentroid = line(handles.time_axis, handles.TumourXcentroid);
TumourYcentroid = line(handles.time_axis, handles.TumourYcentroid);
DiaphragmYbound = line(handles.time_axis, handles.DiaphragmYbound);

% Adjust line properties (functional)
set(TumourXcentroid                          , ...
  'Color'           , [0 0 .7]    );
set(TumourYcentroid                            , ...
  'Marker'          , '.'         , ...
  'Color'           , [0 .7 0]  );
set(DiaphragmYbound                            , ...
    'LineStyle'       , '--'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.7 0 0]  );

% Adjust line properties (aesthetics)
set(TumourXcentroid                          , ...
  'LineWidth'       , 2           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [0 0 .7]   );
set(TumourYcentroid                            , ...
  'LineWidth'       , 2           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [0 .7 0]   );
set(DiaphragmYbound                          , ...
  'LineWidth'       , 2           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 0 0]   );


% Add labels
hTitle  = title ('Pixel location of Ctumour and DiaphragmY');
hXLabel = xlabel('Time(s)'       );
hYLabel = ylabel('Displacement(mm)'    );

% Add legend
hLegend = legend([TumourXcentroid, TumourYcentroid, DiaphragmYbound], ...
  'Ctumour-x AP/RL' , ...
  'Ctumour-y CC', ...
  'Diaphragm-y CC');
  
% Adjust font
set( gca                             , 'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel]      ,  'FontName'   , 'AvantGarde');
% set([hLegend, gca]                   , 'FontSize'   , 8           );
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



%% Volume plot of tumour over time
% This plot the volume of tumour, ie, number of pixels with respect to
% time

figure
hold on
% title (' Volume of tumour')
% grid on
% % plot(handles.time_axis, handles.Areatumour, '-r'...
% %     ,handles.time_axis, handles.DiaphragmYbound,'-g')
% plot(handles.time_axis, handles.Areatumour, '-r')
% xlabel('Time (s)')
% ylabel('Volume(Number of pixels)')
% % ylabel('Volume(1.5625mm^2)')
% hold off
Areatumour = line(handles.time_axis, handles.Areatumour);


% Adjust line properties (functional)
set(Areatumour                          , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );

% Adjust line properties (aesthetics)
set(Areatumour                          , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );

% Add labels
hTitle  = title ('Area of tumour (N pixel)');
hXLabel = xlabel('Time(s)'       );
hYLabel = ylabel('Volume(Number of pixels)'    );

% Add legend
hLegend = legend([Areatumour], ...
  'Area of tumour' );
  
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


 
 
 
