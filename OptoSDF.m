function OptoSDF(spikeRasters,preAlignWindow,pulseDur,axisHandle)

for cellNum=1:length(spikeRasters)
    if ~exist('axisHandle','var') || isempty(axisHandle)
        figure('Position',[1092 149 708 761]); hold on
    end
    colormap(parula); cmap=colormap;
    
    %% plot sdf
    conv_sigma=1;
    %     xTickSteps=round(preAlignWindow/50)*10;
    [sdf{1}, ~, rastsem{1}]=conv_raster(spikeRasters{cellNum},conv_sigma);
    % [sdf{2}, ~, rastsem{2}]=conv_raster(spikeRasters{keepChan},conv_sigma);
    hold on;
    
    %plot sem
    box off; %subplot(1,1);hold on; box off;
    patch([1:length(sdf{1}),fliplr(1:length(sdf{1}))],[sdf{1}-rastsem{1},fliplr(sdf{1}+rastsem{1})],...
        cmap(cellNum,:),'EdgeColor','none','FaceAlpha',0.2); % [0.16 0.38 0.27]
    % endAlignPloth=subplot(1);hold on; box off;
    % patch([1:length(sdf{2}),fliplr(1:length(sdf{2}))],[sdf{2}-rastsem{2},fliplr(sdf{2}+rastsem{2})],cmap(22,:),'EdgeColor','none','FaceAlpha',0.1);
    %plot sdfs
    FRploth=plot(gca,sdf{1},'Color',cmap(cellNum,:),'LineWidth',1.8);%[0.16 0.38 0.27]
    
    zeroLoc=preAlignWindow-3*conv_sigma;
    set(gca,'XTick',[zeroLoc-10 zeroLoc zeroLoc+10 zeroLoc+20 zeroLoc+40]); %xTickSteps-(1+3*conv_sigma):xTickSteps:(stop-start-6*conv_sigma));
    set(gca,'XTickLabel',[-10 0 10 20 40]); %-(alignmtt-xTickSteps):xTickSteps:stop-(alignmtt+xTickSteps));
    axis(gca,'tight');
    set(gca,'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    xlabel(gca,'Time (ms)'); %,'FontName','Cambria','FontSize',12);
    ylabel(gca,'Firing rate (spikes/s)'); %,'FontName','Cambria','FontSize',12);
    
    % draw opto stim bar
    currylim=get(gca,'YLim');
    OptoStimh=patch([preAlignWindow-(3*conv_sigma), preAlignWindow-(3*conv_sigma),...
        preAlignWindow-(3*conv_sigma)+pulseDur-1, preAlignWindow-(3*conv_sigma)+pulseDur-1], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
    set(gca,'Color','white','FontSize',18,'FontName','calibri');
    %legend
    legend([FRploth,OptoStimh],{'Average firing rate','Optical stimulation'},'FontSize',12);
    legend('boxoff')
end