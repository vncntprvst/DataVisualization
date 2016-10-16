function [clusterClasses,lineSelecIdx,reclassOption]=InteractiveClassification(waveforms,clusterClasses,viewClasses)

% class -1 is hidden
% class -10 cannot change
% viewClasses=[viewClasses,-10];
setClass=0;
reclassOption=0;
switch nargin
    case 0
        lineH=findobj(gca,'Type', 'line');
        if ~isempty([lineH.YData])
            waveforms=reshape([lineH.YData],size([lineH.YData],2)/size(lineH,1),size(lineH,1));
            waveforms=waveforms';
        elseif ~exist('waveforms','var')
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
visibleLines=find(cellfun(@(x) strcmp(x,'on'), {lineH.Visible}));
% [lineH(~visibleLines).Visible]=deal('on');
%draw selection line / rectangle
lineSelecIdx=SelectLines(waveforms);

if sum(lineSelecIdx)>0
    disp(['crossed waveform(s) ' num2str(find(lineSelecIdx'))]);
    try
        disp(['with tag(s) ' num2str(viewClasses(lineSelecIdx'))]);
    catch
    end
%     disp(['crossed waveform(s) ' num2str(find(lineSelecIdx')) ', Tag(s) ']);
%     try
%     {lineH(flip(lineSelecIdx)).Tag}
% %         figure;plot(waveforms(lineSelecIdx,:)');hold on
% %         lineWF=fliplr(reshape([lineH(visibleLines).YData],...
% %     size([lineH(visibleLines).YData],2)/size(lineH(visibleLines),...
% %     1),size(lineH(visibleLines),1)));
% % lineWF=lineWF';%one waveform per row
% %         plot(lineWF(lineSelecIdx,:)')
%     catch
%     end
%     lineSelecIdx(clusterClasses==-10)=0;
    %show selection
    set(lineH(visibleLines(flip(lineSelecIdx))),'Color',[0.7 0.5 0.2]);
    uistack(lineH(visibleLines(flip(lineSelecIdx))),'top');
    uistack(findobj(gca,'Type', 'patch'),'top');
    %% actions
    %set class
    prompt={'Set waveform class value','Batch reclassify? (0/1)'};
    name='Waveform classification';
    numlines=1;
    defaultanswer={'0','0'};
    answers=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answers)
        setClass=str2double(answers{1});
        reclassOption=str2double(answers{2});
        %         if ~isempty(setClass)
        clusterClasses(lineSelecIdx)=setClass;
        % disappear
        % set(lineH(flip(lineSelecIdx)),'Visible','off')
        set(lineH(visibleLines(flip(~ismember(clusterClasses,viewClasses)))),'Visible','off')
    else
        set(lineH(visibleLines(flip(lineSelecIdx))),'Color',[0 0 0 0.2]);
        %         end
    end
end
% delete
% delete(lineH(flip(lineSelecIdx))) % delete line
% drawnow()