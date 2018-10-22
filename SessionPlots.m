% Recording session plots
cd 'D:\Data\Ephys\vIRt22\vIRt22_1016\SpikeSortingFolder\vIRt22_2018_10_16_20_36_04_5600_50ms1Hz10mW_1_1_export'
fName='vIRt22_2018_10_16_20_36_04_5600_50ms1Hz10mW_1_1_export.result.hdf5';
traces = memmapfile(['../' regexp(fName,'\S+?(?=\.\w+\.\w+$)','match','once') '.dat'],'Format','int16');
spikes=LoadSpikeData(fName,traces);
%% check that templates started at 0 before adding garbage spikes

%% add TTLs

%% add voltage scaling factor and sampling rate

% plots
spikes.unitID=double(spikes.unitID);
% find most frequent units
[unitFreq,unitIDs]=hist(spikes.unitID,unique(spikes.unitID));
[unitFreq,freqIdx]=sort(unitFreq','descend'); unitIDs=unitIDs(freqIdx);

templateNum=64
% Plot ISI
isis = double(diff(spikes.times(spikes.unitID==templateNum)));
figure; hist(isis)

% plot the amplitude
figure
plot(spikes.times(spikes.unitID==templateNum),spikes.amplitude(spikes.unitID==templateNum), '.')

% Plot waveform
figure
plot(mean(spikes.waveforms(spikes.unitID==templateNum,:)))
