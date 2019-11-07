% Import csv file
[fileName,pathName] = uigetfile('*.csv','Select CSV file');
delimiter = ',';
startRow = 2;

formatSpec = '%f%f%[^\n\r]';
fileID = fopen(fullfile(pathName,fileName),'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);

runData = table(dataArray{1:end-1}, 'VariableNames', {'Value','TimestampTimeOfDayTotalMilliseconds'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

% Adjust values
runData.Value=runData.Value-runData.Value(1);
runData.TimestampTimeOfDayTotalMilliseconds=runData.TimestampTimeOfDayTotalMilliseconds-runData.TimestampTimeOfDayTotalMilliseconds(1);

%Plot
folderName=regexp(pathName,['(?<=\' filesep ')\w+' '(?=\' filesep '$)'],'once','match');
figTitle=[folderName ' wheel training'];
figure('name',figTitle);
plot(runData.TimestampTimeOfDayTotalMilliseconds/60000,runData.Value*0.03,...
    'k','Linewidth',2); %3cm off wheel center 
title({figTitle;'Distance run'})
xlabel('Time (mn)')
ylabel('Distance (m)')
