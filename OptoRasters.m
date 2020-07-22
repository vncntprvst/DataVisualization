function OptoRasters(alignedRasters,alignTime,pulseDur,IPI,axisHandle)
%% plots laser stim aligned rasters
% Input: alignedRasters - cell array with stim aligned rasters
%         alignTime - time used to align responses
%         pulseDur - pulse duration 
%% single neuron plots
for cellNum=1:length(alignedRasters)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761])
    end
%     colormap(parula); % cmap=colormap;
    colormap(flipud(gray));
    %% single neuron raster
    imagesc(logical(alignedRasters{cellNum})); %
    % imagesc(MeanChan);
    xlabel('Time (ms)');
    ylabel('Stimulation#','FontSize',12); %'FontWeight','bold'
    % draw alignment bar
    currylim=get(gca,'YLim');
    %     currxlim=get(gca,'XLim');%midl=round(currxlim(2)/20)*10;
    %     set(gca,'XTick',preAlignWindow:50:max(get(gca,'xlim')));
    %     set(gca,'XTickLabel',0:50:max(get(gca,'xlim'))-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');
%     set(gca,'XLim',[0.5 preAlignWindow+60.5],'XTick',0:10:preAlignWindow+60);
    set(gca,'XTick',0:50:size(alignedRasters{cellNum},2)); %0:100:preAlignWindow+200); %'XLim',[preAlignWindow-250.5 preAlignWindow+250.5]
    set(gca,'XTickLabel',(0:50:size(alignedRasters{cellNum},2))-alignTime,'TickDir','out'); %'FontSize',10,'FontName','calibri'
    
    %% opto stim patch
%     patch([alignTime alignTime alignTime+pulseDur alignTime+pulseDur], ...
%         [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%         [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
%     alignTime=alignTime+IPI-1;
%     patch([alignTime alignTime alignTime+pulseDur alignTime+pulseDur], ...
%         [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%         [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
% 
%     set(gca,'Color','white','FontSize',12,'FontName','Helvetica');
    %     title('Neural response to 10mW stimulation intensity, aligned to stimulation onset');
    % hcb = colorbar('southoutside');
    % hcb.Label.String = 'z-scored firing rate';
    
end

%% plot raster showing all channels
% figure('Position',[1050 120 750 790]);
% colormap bone;
% % subplot(2,1,1)
% meanChan=cellfun(@(x) conv_raster(x,1),spikeRasters,'UniformOutput',false);
% meanChan=cell2mat(meanChan');
% imagesc(zscore(meanChan,[])); %
% % imagesc(MeanChan);
% xlabel('Time (ms)');
% ylabel('Channels','FontWeight','bold','FontSize',12);
% % draw alignment bar
% currylim=get(gca,'YLim');
% currxlim=get(gca,'XLim');midl=round(currxlim(2)/20)*10;
% set(gca,'XTick',[midl-preAlignWindow/2 midl midl+postAlignWindow/2]);
% set(gca,'XTickLabel',[-preAlignWindow/2 0 postAlignWindow/2]);
% %opto stim patch
% patch([midl-1 midl-1 midl+1 midl+1], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.1);
%
% % patch([repmat(midl-3,1) repmat(midl+3,1)], ...
% %     [[0 currylim] fliplr([0 currylim])], ...
% %     [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
% title('Neural response to 10mW stimulation intensity, aligned to stimulation onset');
% hcb = colorbar('southoutside');
% hcb.Label.String = 'z-scored firing rate';

% meanChan=cellfun(@(x) conv_raster(x),spikeRasters,'UniformOutput',false);
% meanChan=cell2mat(meanChan);
% subplot(2,1,2)
% imagesc(zscore(meanChan,[]));
% % imagesc(MeanChan);
% xlabel('Time');
% ylabel('Channel','FontWeight','bold','FontSize',12);
% % draw alignment bar
% currylim=get(gca,'YLim');
% currxlim=get(gca,'XLim');midl=round(currxlim/2);
% set(gca,'XTick',[midl-500 midl midl+500]);
% set(gca,'XTickLabel',[-500 0 500]);
% patch([repmat(midl-3,1,2) repmat(midl+3,1,2)], ...
%     [[0 currylim(2)] fliplr([0 currylim(2)])], ...
%     [0 0 0 0],[0.8 0 0],'EdgeColor','none','FaceAlpha',0.8);
% title('Neural response to 80% stimulation intensity, aligned to stimulation onset');
% hcb = colorbar('southoutside');
% hcb.Label.String = 'z-scored firing rate';

