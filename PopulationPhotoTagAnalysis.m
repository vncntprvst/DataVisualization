% Population Phototagging analysis
clearvars
subjectDir=cd;

% List session directories
expDirs= dir([subjectDir filesep '*' filesep 'Analysis']);
expDirs = expDirs([expDirs.isdir]);
expDirs = expDirs(~cellfun(@(folderName) any(strcmp(folderName,...
    {'.','..'})),{expDirs.name}));

for sessionNum=1:numel(expDirs)
    dataDir=fullfile(expDirs(sessionNum).folder,expDirs(sessionNum).name);
    try
        load(fullfile(dataDir,[expDirs(sessionNum).name '_ephys.mat']));
        load(fullfile(dataDir,[expDirs(sessionNum).name '_pulses.mat']));
        if isempty(pulses.TTLTimes)
            load(fullfile(dataDir,[expDirs(sessionNum).name '_export_trial.mat']));
            pulses.TTLTimes=times;
            clearvars times
        end
        %load traces
        traces = memmapfile(fullfile(dataDir,[expDirs(sessionNum).name '_traces.bin']),'Format','single');
        traces = traces.Data;
        ephys.traces=reshape(traces,[ephys.recInfo.numRecChan ephys.recInfo.dataPoints]);
        clearvars traces
    catch
        disp(['error loading data from ' dataDir])
        continue
    end
    handpick=[]; %46 %vIRt47_0803_5900 : 54
    ephys = SelectUnits(ephys,'quality',handpick); %quality frequency
    ephys = OrderByDepth(ephys);
    
    if ~isfield(pulses,'duration'); pulses.duration=0.01; end
    taggedCells = FindPhototagged(ephys,pulses);
    ephys.selectedUnits  = find(taggedCells);
    PhotoTagPlots(ephys,pulses); 
       
end





