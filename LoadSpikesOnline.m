function spikes = LoadSpikesOnline(source)
%LOADSPIKESONLINE Load an online-sorted (REX) session as a spikes-only model.
%   spikes = LOADSPIKESONLINE(matFile) reads a processed REX session file
%   (e.g. D:\CB paper\0-data\SLData\processed\<animal>\<tag>_REX.mat) that
%   stores ONLY online spike times (no continuous voltage trace and no
%   waveforms) and returns a model the SpikeVisualizationApp can display.
%
%   The app shows the ISI histogram and autocorrelogram from the spike times;
%   the waveform and raw-trace panels are hidden automatically because this
%   source carries no waveforms / data file.
%
%   REX layout used (all times in milliseconds, one row per trial):
%       allspk     [nTrials x maxSpikes]  within-trial spike times (NaN-padded)
%       allspklen  [1 x nTrials]          spike count per trial
%       allstart   [1 x nTrials]          absolute trial-start time (ms)
%       allspkchan [1 x nTrials]          online spike channel (the "cluster")
%   Absolute spike time (ms) = allstart(tr) + allspk(tr, 1:allspklen(tr)).
%
%   See also LoadSpikesPhy, SpikeVisualizationApp.

    arguments
        source (1,1) string {mustBeFile}
    end

    S = load(source);
    required = ["allspk", "allstart", "allspklen"];
    if ~all(isfield(S, required))
        error("LoadSpikesOnline:badFile", ...
            "%s is not a REX session file (needs allspk/allstart/allspklen).", ...
            source);
    end

    nTr = size(S.allspk, 1);
    maxSpk = size(S.allspk, 2);
    if isfield(S, "allspkchan") && numel(S.allspkchan) == nTr
        chanPerTrial = double(S.allspkchan(:));
    else
        chanPerTrial = ones(nTr, 1);
    end

    times = cell(nTr, 1);
    clus = cell(nTr, 1);
    for tr = 1:nTr
        n = S.allspklen(tr);
        if ~isfinite(n) || n <= 0
            continue
        end
        n = min(round(n), maxSpk);
        t = S.allspk(tr, 1:n);
        t = t(~isnan(t));
        if isempty(t)
            continue
        end
        abs_ms = double(S.allstart(tr)) + double(t(:));   % absolute ms
        times{tr} = abs_ms;
        c = chanPerTrial(tr);
        if ~isfinite(c)
            c = 0;
        end
        clus{tr} = repmat(c, numel(abs_ms), 1);
    end
    times = vertcat(times{:});
    clus = vertcat(clus{:});
    if isempty(times)
        error("LoadSpikesOnline:noSpikes", "%s contains no spikes.", source);
    end

    % Times are already in ms -> treat as a 1 kHz sample clock so downstream
    % ISI/ACG math (which divides by samplingRate) yields milliseconds.
    fs = 1000;
    [times, order] = sort(times);
    clus = clus(order);

    [~, base] = fileparts(source);
    spikes.sourceType = "online";
    spikes.name = string(base);
    spikes.sourcePath = string(source);
    spikes.samplingRate = fs;
    spikes.spikeTimes = round(times);
    spikes.clusters = clus;
    spikes.clusterIds = unique(clus);
    spikes.numSamples = max(spikes.spikeTimes);
    spikes.numChannels = 1;
    spikes.dataChannel = 0;
    spikes.gainToUV = 1;                 % no voltage / no waveforms
    spikes.templates = [];
    spikes.amplitudes = [];
    spikes.waveforms = single.empty(numel(times), 0);
    spikes.waveformWindow = [40 40];
    spikes.pcFeatures = single.empty(numel(times), 0, 0);
    spikes.pcFeatureInd = [];
    spikes.clusterTable = table(spikes.clusterIds(:), VariableNames="cluster_id");
end
