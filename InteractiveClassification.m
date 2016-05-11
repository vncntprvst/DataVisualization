function InteractiveClassification(waveforms,clusterClasses,viewClasses)

switch nargin
    case 0
        if ~exist('waveforms','var')
            % example waveforms, one per row
            waveforms=[(sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11));...
                (sin(-10:0)+rand(1,11))];
            
            % plot WFs
            figure;
            plot(waveforms');
        end
    case 1
        clusterClasses=ones(size(waveforms,1),1);
        viewClasses=unique(clusterClasses);
end

%get line handles
lineH=findobj(gca,'Type', 'line');%in reverse order

%draw selection line / rectangle
lineSelecIdx=SelectLines(waveforms);

%% actions
%set class 
prompt='Set waveform class value';
name='Waveform classification';
numlines=1;
defaultanswer={'0'};
setClass=str2double(inputdlg(prompt,name,numlines,defaultanswer));
clusterClasses(lineSelecIdx)=setClass;

% disappear
% set(lineH(flip(lineSelecIdx)),'Visible','off')
set(lineH(flip(~ismember(clusterClasses,viewClasses))),'Visible','off')

% delete
% delete(lineH(flip(lineSelecIdx))) % delete line
% drawnow()