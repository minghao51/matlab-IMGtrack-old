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
handles.TumourXcentroid = handles.CtumourX*handles.info{1}.PixelSpacing(1);
handles.TumourXcentroid = handles.TumourXcentroid - nanmean(handles.CtumourX)*handles.info{1}.PixelSpacing(1); %Mean ignoring NaN values
handles.TumourYcentroid = handles.CtumourY*handles.info{1}.PixelSpacing(2);
handles.TumourYcentroid = handles.TumourYcentroid - nanmean(handles.CtumourY)*handles.info{1}.PixelSpacing(1);
handles.DiaphragmYbound = handles.DiaphragmY*handles.info{1}.PixelSpacing(2);
handles.DiaphragmYbound = handles.DiaphragmYbound - nanmean(handles.DiaphragmY)*handles.info{1}.PixelSpacing(1);

%% Plot of Centroid of tumour and Diaphragm coordinates
% This plot the appropriate data, with the respectice time in s and
% displacement in mm

figure
hold on
% title(' pixel location of Ctumour and diaphragm');
% grid on
% set(gca,'GridLineStyle','-');
% grid minor;
% plot(handles.time_axis,handles.TumourXcentroid,'-r'...
%     ,handles.time_axis,handles.TumourYcentroid,'-b'...
%     ,handles.time_axis,handles.DiaphragmYbound,'-g')
% xlabel('Time (s)')
% ylabel('Displacement(mm)')
% 
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
hTitle  = title ('Pixel location of Ctumour and diaphragm');
hXLabel = xlabel('Time(s)'       );
hYLabel = ylabel('Displacement(mm) [inferior\caudal]'    );

% Add legend
hLegend = legend([TumourXcentroid, TumourYcentroid, DiaphragmYbound], ...
  'Ctumour-x AP/RL' , ...
  'Ctumour-y CC', ...
  'Diaphragm-y CC');
  
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

%% Motion Speed
handles.TumourDisplacement = sqrt(handles.TumourXcentroid.^2 ...
    + handles.TumourYcentroid.^2);
handles.Speed = diff(handles.TumourDisplacement);
handles.time_axis_speed = ((0+temporaltime)/2:temporaltime:(endpoint));
figure,plot (handles.time_axis_speed,handles.Speed);

%% Volume plot of tumour over time
% This plot the volume of tumour, ie, number of pixels with respect to
% time.
% The Expected volume of the phantom bottle is 1607.3 +15 pixel
EVolume = 1607.3;
EVolumeUncertainty = 150;

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
Expected = line([0, 30],[EVolume, EVolume] );
ExpectedL = line([0, 30],[EVolume-EVolumeUncertainty, EVolume-EVolumeUncertainty]  );
ExpectedU = line([0, 30], [EVolume+EVolumeUncertainty, EVolume+EVolumeUncertainty]);

% Adjust line properties (functional)
set(Areatumour                          , ...
  'Marker'          , '.'         , ...
  'Color'           , [.5 .5 .5]  );
set(Expected                        , ...
  'LineStyle'       , '--'        , ...
  'Color'           , 'r'         );
set(ExpectedL                        , ...
  'LineStyle'       , '-.'        , ...
  'Color'           , [0 .5 0]    );
set(ExpectedU                        , ...
  'LineStyle'       , '-.'        , ...
  'Color'           , [0 .5 0]    );

% Set the axis limits
% axis([0 30 1570 1600]);

% Adjust line properties (aesthetics)
set(Areatumour                          , ...
  'LineWidth'       , 1           , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 6           , ...
  'MarkerEdgeColor' , [.2 .2 .2]  , ...
  'MarkerFaceColor' , [.7 .7 .7]  );
set(Expected                        , ...
  'LineWidth'       , 1.5         );
set(ExpectedL                       , ...
  'LineWidth'       , 1.5         );
set(ExpectedU                     , ...
  'LineWidth'       , 1.5         );

% Add labels
hTitle  = title ('Area of tumour (N pixel)');
hXLabel = xlabel('Time(s)'       );
hYLabel = ylabel('Volume(Number of pixels)'    );

% Add legend
hLegend = legend([Areatumour, Expected, ExpectedL], ...
  'Area of tumour' ,...
  'Measured Area',...
  'Uncertainty');
  
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

