
[Behavior,Performance]=processBehaviorData;
keepSessions=find(mean(vertcat(Performance.CorrectSides),2)>5);
Performance=Performance(keepSessions);
Behavior=Behavior(keepSessions);
% [~,orderDate]=sort([Behavior(keepSessions).fileRecordingDate]);
[~,orderDate]=sort(([Behavior.fileRecordingDate]-...
    datetime(year(min([Behavior.fileRecordingDate])),1,1))+...
    seconds([Behavior.fileStartTime_ms]/1000));
Performance=Performance(orderDate);
Behavior=Behavior(orderDate);
% average side bias
avgpcCorrectSides=vertcat(Performance.CorrectSides)./repmat(sum(vertcat(Performance.CorrectSides),2),1,2)*100;

% Evolution of performance
evolPerf=[Performance.overallPerf];
%calcute bootstrapped 95% confidence intervals

% Get d' and decision criterion (~bias)
[dprime,crit] = SigDetecPerformance([Performance.hitRate],[Performance.falseAlarm]);

%find trial regime (randomized or block). If blocks, add patch on block
%session, and add block size on graph
sessions.BlockSize=nan(max(size(keepSessions)),1);
for sessionNum=1:max(size(keepSessions))
    blockSize=unique(diff(find(diff(Behavior(sessionNum).trials.trialType)~=0)));
    if max(size(blockSize))==1
        sessions.BlockSize(sessionNum)=blockSize;
    else
        sessions.BlockSize(sessionNum)=0; %random
    end
end
sessions.types.blocks=bwlabel(sessions.BlockSize);
sessions.types.blocksIdx=find(bwlabel(sessions.BlockSize));
sessions.types.random=bwlabel(sessions.BlockSize==0);
sessions.types.randomIdx=find(bwlabel(sessions.BlockSize==0));

%plots
figure('Name','','NumberTitle','off','position',[1000 215 800 750])
colormap lines;
cmap = colormap(gcf);

% plot side bias
subplot(2,2,1); hold on;
try %if violin plot is available
    % rather than a square plot, make it thinner
    violinPlot(avgpcCorrectSides(:, 1), 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 0, ...
        'color',  mat2cell(cmap(1, : ), 1));
    
    violinPlot(avgpcCorrectSides(:, 2), 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 0, ...
        'color',  mat2cell(cmap(2, : ), 1));
    set(gca, 'xtick', [0.6 1.4], 'xticklabel', {'Left','Right'}, 'xlim', [0.2 1.8]);
    
    % add significance stars for each bar (Not needed here)
        xticks = get(gca, 'xtick');
    %     for side = 1:2,
    %         [~, pval] = ttest(avgpcCorrectSides(:, side));
    % %         yval = max(avgpcCorrectSides(:, side)) * 1.2; % plot this on top of the bar
    %         yval = 6; % plot below
    %         mysigstar(gca, xticks(side), yval, pval);
    %         % if mysigstar gets just 1 xpos input, it will only plot stars
    %     end
    
    % significance star for the difference
    [~, pval] = ttest(avgpcCorrectSides(:, 1), avgpcCorrectSides(:, 2));
    % if mysigstar gets 2 xpos inputs, it will draw a line between them and the
    % sigstars on top
    mysigstar(gca, xticks, max(max(avgpcCorrectSides)) * 1.1, pval);
    
catch
    bar(mean(avgpcCorrectSides));
    hold on
    errorbar(mean(avgpcCorrectSides),std(avgpcCorrectSides,1),'kx','LineWidth',1)
    set(gca,'xticklabel',{'Left','Right'})
    set(gca,'Color','white','TickDir','out')
end
ylabel('Percentage Side Choice')
title('Correct answers, Left vs Right')
    
% Success rate
subplot(2,2,2)
plot(evolPerf,'LineWidth',1.5,'Color',cmap(4,:))
hold on
plot(1:max(size(evolPerf)),0.75*ones(1,max(size(evolPerf))),'k--','LineWidth',1.5)
axis(gca,'tight'); box off;
set(gca,'Color','white','TickDir','out','ylim',[0 1])
xlabel('Session')
ylabel('Evolution of performance level')
title('Success rate')

% Performance
subplot(2,2,3:4); hold on;
for sessTypeNum=1:max(sessions.types.random)
plot(find(sessions.types.random==sessTypeNum),dprime(sessions.types.random==sessTypeNum),'LineWidth',1.5,'Color',cmap(5,:),'Marker','o','MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.9 0.5 0.1],'MarkerSize',10)
end
for sessTypeNum=1:max(sessions.types.blocks)
plot(find(sessions.types.blocks==sessTypeNum),dprime(sessions.types.blocks==sessTypeNum),'LineWidth',1.5,'Color',cmap(5,:),'Marker','o','MarkerEdgeColor','k',...
                'MarkerFaceColor',[0.1 0.5 0.9],'MarkerSize',10)
end
plot(1:max(size(dprime)),1.5*ones(1,max(size(dprime))),'k--','LineWidth',1.5)
axis(gca,'tight'); box off;
set(gca,'Color','white','TickDir','out')
set(gca,'ylim',[floor(min([0 get(gca,'ylim')])) max([2 get(gca,'ylim')])]);
try
set(gca,'xtick',linspace(1,size(Behavior,2),size(Behavior,2)));
set(gca,'xticklabel',cellfun(@(recdate) datestr(recdate,'mmm-dd'), {Behavior.fileRecordingDate},'UniformOutput',false))
catch
end
xlabel('Session')
ylabel('Performance (d'')')
if ~isempty(sessTypeNum)
    legend({'Random trials','Block trials','Discrimination performance criterion'},'location','SouthEast');
else
    legend({'Random trials','Discrimination performance criterion'},'location','SouthEast');
end

title('Texture detection performance')

