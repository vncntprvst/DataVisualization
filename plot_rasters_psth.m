function plot_rasters_psth(varargin)
figure('Position',[1050 120 750 790]);
cmap=colormap;
for dataRow=1:size(varargin{:},1)
    args=varargin{:}{dataRow};
    if size(args,2)==1
        % Figure options is a structure that specifies limits, alignments, style, legends, etc
        % alignment specifications is a cell with 4 parameters:
        %         interval pre-alignment
        %         interval post-alignment
        %         alignment bar width
        %         alignment bar color
        % figure style is a cell with 3 parameters:
        %         colormap
        %         x tick spacing
        %         x tick labels
        % legends is a cell with 3 parameters:
        %         x axis label
        %         y axis label
        %         title
        
        figureOptions=struct('alignSpecs',...
            {size(rasterData,2)/2+1;size(rasterData,2)/2+1;1;'white'},...
            'figureStyle',{'parula';1:10:size(rasterData,2)+1;-(size(rasterData,2)/2+1)/1000:(size(rasterData,2)/2+1)/1000},...
            'legends',{'Time (s)'; 'Trials'; {['Neuron # ' num2str(size(rasterData(1)))],...
            'response aligned to midpoint'}});
    else
        figureOptions=args{2};
    end
    
    % find rasterplot type
    if iscell(args{1}) && numel(args{1})==2
        plotType='both';
        signalData=args{1}{1};
        spikeData=args{1}{2};
    elseif iscell(args{1})
        plotType='spikes';
        spikeData=args{1};
    else
        plotType='traces';
        signalData=args{1};
    end
    
    subplot(4,size(varargin{:},1),dataRow:size(varargin{:},1):dataRow+2*(size(varargin{:},1)))
    colormap(figureOptions.figureStyle{1});
    
    if contains(plotType,'traces') | contains(plotType,'both')
        %# method 1
        imagesc(signalData); colorbar('northoutside')%     imagesc(zscore(spikeRasters,[])); % spikeArray
        hold on
        % imagesc(MeanChan);
    end
    if contains(plotType,'spikes') | contains(plotType,'both')
        %#method 2:  plot spike indices
        if iscell(spikeData) % list of time indices
            trialWithSpikes=find(~cellfun('isempty', spikeData));
            numSpikes=cellfun(@(x) size(x,1), spikeData(trialWithSpikes));
            spikeTimes=cellfun(@(x) x(1:2:end), spikeData(trialWithSpikes),'UniformOutput',false);
            try %for calcium events
                spikePeak=cellfun(@(x) x(2:2:end), spikeData(trialWithSpikes),'UniformOutput',false);
            catch
                spikePeak=[];
            end
            [indx,indxPeak,indy]=deal(NaN(sum(numSpikes),1));
            spikeNum=1;
            for trialNum=1:numel(trialWithSpikes)
                for thatTrialSpikeNum=1:numSpikes(trialNum)
                    indy(spikeNum)=trialWithSpikes(trialNum);
                    indx(spikeNum)=spikeTimes{trialNum}(thatTrialSpikeNum);
                    if ~isempty(spikePeak)
                        indxPeak(spikeNum)=spikePeak{trialNum}(thatTrialSpikeNum);
                    end
                    spikeNum=spikeNum+1;
                end
            end
        elseif ismatrix(spikeData)
            [indy, indx] = ind2sub(size(spikeData),find(spikeData)); %find row and column coordinates of spikes if it's
        end
        if ~isnan(sum(indxPeak)) %extended rasters for calcium events
            % simple plot method
            %         plot([indx';indxPeak'],[indy';indy'],'color','k','LineStyle','-','LineWidth',2);
            % patch method
            patch('Faces',reshape(1:numel(indx)*4,[4,numel(indx)])',...
                'Vertices',reshape([indx, indx, indxPeak, indxPeak;...
                indy-0.5, indy+0.5, indy+0.5, indy-0.5]',4*numel(indx), []),...
                'FaceVertexCData',repmat([0;0;6;6],numel(indx),1),...
                'FaceColor','interp','EdgeColor','none','FaceAlpha',0.8);  % 'FaceColor',[0.4 0.1 0.3]
        else
            plot([indx';indx'],[indy'-0.5;indy'+0.5],'color','k','LineStyle','-','LineWidth',2); % plot rasters
        end
        %     set(gca,'Ydir','reverse')
    end
    
    ylabel(figureOptions.legends{2},'FontWeight','bold','FontSize',12);
    
    currylim=get(gca,'YLim');
    %     currxlim=get(gca,'XLim');
    set(gca,'XTick',figureOptions.figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figureOptions.figureStyle{3});
    
    % draw alignment bar
    patch([repmat(figureOptions.alignSpecs{1},1,2) repmat(figureOptions.alignSpecs{1}+figureOptions.alignSpecs{3},1,2)], ...
        [[0 currylim(2)] fliplr([0 currylim(2)])], ...
        [0 0 0 0],figureOptions.alignSpecs{4},'EdgeColor','none','FaceAlpha',0.3);
    box off; % axis tight
    title(figureOptions.legends{3});
    %     hcb = colorbar('eastoutside');
    %     hcb.Label.String = figureOptions.legends{4};
    
    
    %% plot psth
    subplot(4,size(varargin{:},1),dataRow+3*(size(varargin{:},1)))
    %     for bar plots:
    %     barPlot=bar(mean(signalData));
    %     %     barPlot.FaceColor=[0.1 0.4 0.8];
    %     barPlot.EdgeColor='none';
    %     barPlot.BarWidth=1;
    
    % for line plots:
    plot(mean(signalData));
    set(gca,'ylim',[-0.4 0.4]);
    currylim=get(gca,'YLim');
    % draw sem
    patch([1:size(signalData,2) flip(1:size(signalData,2))], ...
        [mean(signalData)+(std(signalData)/sqrt(size(signalData,1)))...
        fliplr(mean(signalData)-(std(signalData)/sqrt(size(signalData,1))))], ...
        cmap(1,:),'EdgeColor','none','FaceAlpha',0.3);
    % draw alignment bar
    patch([repmat(figureOptions.alignSpecs{1},1,2)-0.2 ...
        repmat(figureOptions.alignSpecs{1}+0.2,1,2)], ...
        [currylim fliplr(currylim)], ...
        'k','EdgeColor','none','FaceAlpha',0.5);
    set(gca,'XTick',figureOptions.figureStyle{2},'TickDir','out');
    set(gca,'XTickLabel',figureOptions.figureStyle{3});
    box off; axis tight
    %     set(gca,'ylim',[0 0.5]);
    ylabel('Mean z-scored \DeltaF/F0','FontWeight','bold','FontSize',12,'Interpreter','tex');
    xlabel(figureOptions.legends{1},'FontWeight','bold','FontSize',12);
end


