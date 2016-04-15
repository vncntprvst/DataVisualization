
[Behavior,Performance]=processBehaviorData;
keepSessions=find(mean(vertcat(Performance.ReportSides),2)>30);
[~,orderDate]=sort([Behavior(keepSessions).fileRecordingDate]);
Performance=Performance(keepSessions(orderDate));
% average side bias
avgpcReportSides=vertcat(Performance.ReportSides)./repmat(sum(vertcat(Performance.ReportSides),2),1,2)*100;

% Evolution of performance
evolPerf=[Performance.overallPerf];
%calcute bootstrapped 95% confidence intervals

%find trial regime (randomized or block). If blocks, add patch on block
%session, and add block size on graph

% plot side bias

figure('Name','','NumberTitle','off','position',[1000 215 800 750])
colormap lines;
cmap = colormap(gcf);

subplot(1,2,1)
bar(mean(avgpcReportSides));
hold on
errorbar(mean(avgpcReportSides),std(avgpcReportSides,1),'kx','LineWidth',1)
set(gca,'xticklabel',{'Left','Right'})
set(gca,'Color','white','TickDir','out')
ylabel('Percentage Side Choice')
title('Correct answers, Left vs Right')

subplot(1,2,2)
plot(evolPerf,'LineWidth',1.5,'Color',cmap(4,:))
axis(gca,'tight'); box off;
set(gca,'Color','white','TickDir','out','ylim',[0 1])
xlabel('Time (mn)')
ylabel('Evolution of performance level')
title('Success rate')


