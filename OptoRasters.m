function OptoRasters(spikeRasters,preAlignWindow,pulseDur,axisHandle)

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

%% single neuron plots
for cellNum=1:length(spikeRasters)
    %find(mean(meanChan)==max(mean(meanChan)));
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761])
    end
    colormap(parula); % cmap=colormap;
    
    %% single neuron raster 
    imagesc(spikeRasters{cellNum}); %
    
    % imagesc(MeanChan);
    xlabel('Time (ms)');
    ylabel('Stimulation#','FontSize',12); %'FontWeight','bold'
    % draw alignment bar
    currylim=get(gca,'YLim');
%     currxlim=get(gca,'XLim');%midl=round(currxlim(2)/20)*10;
    set(gca,'XTick',preAlignWindow:50:max(get(gca,'xlim')));
    set(gca,'XTickLabel',0:50:max(get(gca,'xlim'))-preAlignWindow,'FontSize',10,'FontName','calibri','TickDir','out');
    
    %opto stim patch
    patch([preAlignWindow preAlignWindow preAlignWindow+pulseDur-1 preAlignWindow+pulseDur-1], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    set(gca,'Color','white','FontSize',18,'FontName','calibri');
%     title('Neural response to 10mW stimulation intensity, aligned to stimulation onset');
    % hcb = colorbar('southoutside');
    % hcb.Label.String = 'z-scored firing rate';
  
end

