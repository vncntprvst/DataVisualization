cd('E:\Data\export\PrV88_119_OEph_16Ch');
% fileName='PrV88_119AN_OptoStim5_OEph_16Ch';
fileName='PrV88_119AN_WakingOptoStim_OEph_16Ch';
stimParams='2mW 2ms 5Hz';

% LFPdata=rawData;
LFP_Data=importdata([fileName '_LFP_raw.mat']);
trialData=importdata([fileName '_CAR_trials.mat']);

preproc=1;

% set parameters
spectrumParams.Fs=1000;         % sampling frequency
spectrumParams.fpass=[0 100];   % frequency band of interest
spectrumParams.tapers=[5 9];   % taper parameters [TW K].
% TW is the  time-bandwidth product
% K is the number of tapers to be used (less than or equal to 2TW-1).
spectrumParams.pad=2;           % padding factor for fft
% -1 no padding,
% 0 padding to the next highest power of 2

% format data
LFP_Data=single(LFP_Data);
if size(LFP_Data,2)>size(LFP_Data,1)
    LFP_Data=LFP_Data';
end

%% preprocessing - if not done already
if preproc
    % Remove DC offsets and slowly changing components with locdetrend
    % function, using 1s moving window.
    movingWin=[.1 .05]; % [1 0.1]
    LFP_Data=locdetrend(LFP_Data,spectrumParams.Fs,movingWin);
    %     [S,f]=mtspectrumc(LFPdata,params);
    % 	figure; plot_vector(S,f);
    
    % remove 60Hz line noise
    LFP_Data=rmlinesc(LFP_Data,spectrumParams);
    %     [S,f]=mtspectrumc(LFPdata,params);
    % 	figure; plot_vector(S,f);
end

%% Average power spectral density estimation
% first filter the raw LFP trace using a 300-Hz low-pass filter
% (fourth order Butterworth; filter and filtfilt routines)
% and downsample to 663.5-Hz sampling frequency.
% A Welch spectral estimator is then applied to obtain the PSD
% (MATLAB Signal Processing Toolbox, spectrum.welch,
% using a single Hamming window taper and 50% overlapping 200-ms time windows).

spectrumParams.err= 0;          % error calculation
                        % [1 p] - Theoretical error bars;
                        % [2 p] - Jackknife error bars eg [2 0.05];
                        % [0 p] or 0 - no error bars
spectrumParams.trialave=0;      % trial (or channel) average
                        % 1 average over trialData/channels
                        % 0 don't average
movingWin=[0.5 0.05];   % moving window
dBrange=[0 35];
% plot one channel / trial
% [spectrum,timeVals,freqVals] = mtspecgramc(LFPdata(:,1), movingWin, params);
% figure;colormap(jet);
% imagesc(timeVals,freqVals,10*log10(spectrum)'); axis xy;

% each channel / trial separately
[spectrum,timeVals,freqVals] = mtspecgramc(LFP_Data, movingWin, spectrumParams);

figure;colormap(jet);
for chNum=1:size(LFP_Data,2)
    subplot(ceil(size(LFP_Data,2)/4),4,chNum)
    imagesc(timeVals,freqVals,10*log10(spectrum(:,:,chNum))',dBrange); axis xy;
    title(['Channel ' num2str(chNum)])
end
cbH=colorbar;
set(cbH, 'Position', [.92 .11 .03 .8150]); ylabel(cbH, 'Spectrum dB')
invAxH=axes('position',[.1 .1 .85 .85],'visible','off');
xlabel('Time (s)','visible','on','FontName','Helvetica','FontSize',14); 
ylabel('Frequency (Hz)','visible','on','FontName','Helvetica','FontSize',14);
title(['Spectrograms across channels. Antidromic VPM->PrV. ' stimParams],'visible','on','FontName','Helvetica','FontSize',15);
fNameText=text(min(get(gca,'xlim'))-0.1, max(get(gca,'ylim'))+0.04, fileName);
set(fNameText,'FontName','Helvetica','FontSize',6,'Interpreter','none');

% average LFP across channels / trialData
spectrumParams.trialave=1;
movingWin=[0.1 0.01];   % moving window
[spectrum,timeVals,freqVals] = mtspecgramc(LFP_Data, movingWin, spectrumParams);
figure;colormap(jet);
% subplot(10,1,2:10);
imagesc(timeVals*1000,freqVals,10*log10(spectrum)',dBrange); axis xy; cbH=colorbar; ylabel(cbH, 'Spectrum dB');
xlabel('Time (s)','FontName','Helvetica','FontSize',14); 
ylabel('Frequency (Hz)','FontName','Helvetica','FontSize',14);
box off; %axis(gca,'tight'); 
set(gca,'Color','white','TickDir','out','FontName','Helvetica','FontSize',14);
spectrumPosition=get(gca,'position');
% xLims=get(gca,'xlim');
yLims=get(gca,'ylim');
set(gca,'ylim',[-4 yLims(2)]);
fNameText=text(min(get(gca,'xlim'))-10, max(get(gca,'ylim'))+7, fileName);
set(fNameText,'FontName','Helvetica','FontSize',6,'Interpreter','none');
% plot pulse
% subplot(10,1,1);
% axes('position',spectrumPosition,'visible','off')
% set(gca,'position',[pulsePosition(1), pulsePosition(2) spectrumPosition(3) pulsePosition(4)]);
trialData.start(:,3)=(trialData.start(:,2)-(trialData.startClockTime/30));
trialData.end(:,3)=(trialData.end(:,2)-(trialData.startClockTime/30));
% too small!
for TTLNum=1:size(trialData.start,1)
    patch([trialData.start(TTLNum,3):trialData.end(TTLNum,3),...
        fliplr(trialData.start(TTLNum,3):trialData.end(TTLNum,3))],...
        [ones(1,trialData.end(TTLNum,3)-trialData.start(TTLNum,3)+1)*(yLims(1)-1),...
        ones(1,trialData.end(TTLNum,3)-trialData.start(TTLNum,3)+1)*(yLims(1)-2)],...
        [0 0 0],'EdgeColor','black','FaceAlpha',1); %[0.3 0.75 0.93]
end
patch([trialData.start(1,3):trialData.end(end,3),...
    fliplr(trialData.start(1,3):trialData.end(end,3))],...
    [ones(1,trialData.end(end,3)-trialData.start(1,3)+1)*(yLims(1)-2),...
    ones(1,trialData.end(end,3)-trialData.start(1,3)+1)*(yLims(1)-4)],...
    [0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.5);
% set(gca,'xlim',xLims,'ylim',[-4 yLims(2)]);
title(['Average LFP across channels. Antidromic VPM->PrV. ' stimParams],'FontName','Helvetica','FontSize',15);
box off;
% set(gca,'Color','white','XColor','white', 'YColor','white','XTick',[],'YTick',[]);
% set(gca,'Color','white', 'XTick',[],'YTick',[],'Visible','off');
set(gca,'XTickLabel', num2str(cellfun(@(label) str2double(label)*10,get(gca,'XTickLabel'))))
% set(gca,'XTickLabel','')
% Pulse-triggered average
trialData.start=single(trialData.start)/1000;
trialData.end=single(trialData.end)/1000;
trialData.startClockTime=single(trialData.startClockTime);
dBrange=[10 30];
spectrumParams.trialave=1; %Average across trials
win=[1 1]; %1 sec before to 1 sec after cue
movingWin=[0.05 0.005];   % moving window
[spectrum,timeVals,freqVals] = mtspecgramtrigc(LFP_Data(:,1),trialData.start(:,3),win,movingWin,spectrumParams);
figure
imagesc(timeVals-win(1),freqVals,10*log10(spectrum'),dBrange)
axis xy %Flip Y-axis
set(gca,'Color','white','TickDir','out','FontName','Helvetica','FontSize',14);
xlabel('Time (s)','FontName','Helvetica','FontSize',14); 
ylabel('Frequency (Hz)','FontName','Helvetica','FontSize',14);
box off; %axis(gca,'tight'); 
title('Pulse-triggered LFP average','FontName','Helvetica','FontSize',15);
set(gca,'XTick',-1:0.5:1)
cbH=colorbar; ylabel(cbH, 'Spectrum dB');
fNameText=text(min(get(gca,'xlim'))-0.3, max(get(gca,'ylim'))+8, fileName);
set(fNameText,'FontName','Helvetica','FontSize',6,'Interpreter','none');

%% Gamma power estimation
% To estimate changes in gamma power on behavioral time scales,
% a spectrogram is constructed for each LFP trace using the multitaper method
% which estimates the spectral power S(f) in a finite,
% sliding time window by averaging over Fourier transforms
% (discrete Fourier transforms evaluated using the Fast Fourier Transform algorithm on zero-padded data)
% obtained from each of a set of K orthogonal tapers applied to the data.
% The tapers are the first K functions that optimize spectral concentration
% (the tradeoff between broadband and narrowband bias inherent in spectral estimation),
% known as Slepians. We use the Chronux mtspecgramc function, with the following parameters:
% window size, 0.5 s; time step, 50 ms; five tapers; bandwidth 6 Hz.
% From the resulting spectrogram, power in the 45–55 Hz frequency range is averaged
% to obtain a “gamma-50” time series which could then be analyzed as a function of time and spatial location.
% Average power in the 70–85 Hz range is used for “gamma-80”.

%% Spike phase estimation
% To obtain the mean phase angle of spiking relative to ongoing gamma oscillations,
% LFPs are first band-pass filtered (fourth order Chebyshev, r = 0.5,
% MATLAB filter and filtfilt routines; 45–55 Hz for gamma-50, 70–85 Hz for gamma-80)
% before a Hilbert transform is applied to obtain the instantaneous phase angle.
% A histogram of spike counts in each of 10° phase bins is constructed,
% and Rayleigh’s r test used to determine the significance of deviations
% from the uniform distribution (Fisher, 1993 ).




