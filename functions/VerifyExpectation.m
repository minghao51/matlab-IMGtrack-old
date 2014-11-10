function [ handles ] = VerifyExpectation(handles, Tumourpixelmargin, Diahphragmpixelmargin )
%% Verifying Method
% Ensure that recorded data are within expectation
% Check for massive change of disconnection between point based on the
% change of tumour/diaphragm location from one frame to another, and print
% out an error msg.

% Could use mean and variance to explore outlying data as well.

% If Tumourpixelmargin and Diaphragmmargin aren't define, set them as default value
if(exist('Tumourpixelmargin','var')==0), Tumourpixelmargin=10; end
if(exist('Diahphragmpixelmargin','var')==0), Diahphragmpixelmargin=20; end

%Loop through the position of Tumour and compare with one before, the
%chances of position should be within Tumourpixelmargin
% for i = 2: numel(handles.dicomlist)
%     if abs(handles.CtumourX(i-1) - handles.CtumourX(i)) > Tumourpixelmargin
%         fprintf(1, 'Verifying, Exception at %d CtumourX\n',i-1);
%         handles.CtumourX(i-1) = NaN;
%         handles.CtumourY(i-1) = NaN;
%         handles.Areatumour(i,1) = NaN;
%     end
%     if abs(handles.CtumourY(i-1) - handles.CtumourY(i)) > Tumourpixelmargin
%         fprintf(1, 'Verifying, Exception at %d CtumourY\n',i-1);
%         handles.CtumourX(i-1) = NaN;
%         handles.CtumourY(i-1) = NaN;
%     end
%     if handles.DiaphragmY(i-1)-handles.DiaphragmY(i)  > Diahphragmpixelmargin;
%         fprintf(1, 'Verifying, Exception at %d DiaphragmY\n',i-1);
%         handles.DiaphragmY(i) = NaN;
%         handles.CtumourX(i) = NaN;
%         handles.CtumourY(i) = NaN;
%         handles.Areatumour(i,1) = NaN;
%     end
% end


% nanmean(handles.Areatumour)

%Seeing that the data below 1.2 times of the volume of the first frame have more significant trends
thresh = nanmean(handles.Areatumour)*1.4;
%Utilizing arrayfun/cellfun to check for indices of arrays that is below
%threshold.
vidx = arrayfun(@(x) x < thresh, handles.Areatumour);
%Preserving only those values smaller than threshold level
handles.Areatumour = handles.Areatumour(vidx);
handles.CtumourY= handles.CtumourY(vidx);
handles.CtumourX = handles.CtumourX(vidx);
handles.DiaphragmY = handles.DiaphragmY(vidx);
handles.time_axis = handles.time_axis(vidx);

%     % Loop through to check for position changes, not working atm
%     for i = 2: numel(handles.dicomlist)
%         if abs(handles.Areatumour(cnt,1)- handles.Areatumour(cnt-1,1))> handles.Areatumour(cnt,1)
%             handles.CtumourX(i-1) = NaN;
%             handles.CtumourY(i-1) = NaN;
%             handles.Areatumour(i,1) = NaN;
%         end
%     end
