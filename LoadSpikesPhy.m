function spikes = LoadSpikesPhy(phyDir, options)
%LOADSPIKESPHY Load a Phy/KiloSort-format spike sorting into a MATLAB struct.
%   spikes = LOADSPIKESPHY(phyDir) reads the Phy export folder phyDir (the
%   folder that contains params.py, spike_times.npy, spike_clusters.npy, ...)
%   and returns a struct describing every spike, its cluster, the per-spike
%   waveform extracted from the raw data channel, and the sorter PC features.
%
%   spikes = LOADSPIKESPHY(phyDir, Name=Value) accepts:
%       DataChannel  - 0-based index of the real signal channel in data.dat.
%                      Default 0. (These exports store ch0 = real signal and
%                      ch1 = a synthetic duplicate, so ch0 is the default.)
%       PreSamples   - samples kept before each spike time. Default 40.
%       PostSamples  - samples kept after each spike time.  Default 40.
%       ExtractWaveforms - logical, extract per-spike waveforms from data.dat.
%                      Default true. Set false to skip (e.g. data.dat absent).
%
%   Output fields:
%       sourceType     "phy"
%       phyDir         resolved absolute path of the export folder
%       datPath        full path of the raw data file (data.dat)
%       samplingRate   sampling rate in Hz
%       gainToUV       scale factor from raw ADC units to microvolts
%       numChannels    number of channels stored in data.dat
%       numSamples     recording length in samples (per channel)
%       dataChannel    0-based channel used for waveform extraction
%       spikeTimes     nSpikes-by-1 double, spike times in samples
%       clusters       nSpikes-by-1 double, cluster id per spike
%       templates      nSpikes-by-1 double, template id per spike
%       amplitudes     nSpikes-by-1 double, template amplitude per spike
%       waveforms      nSpikes-by-nWfSamples single, ch0 waveform per spike
%       waveformWindow [preSamples postSamples]
%       pcFeatures     nSpikes-by-nPC-by-nLocalChan single, sorter PC features
%       pcFeatureInd   nTemplates-by-nLocalChan, channel index per PC column
%       clusterTable   table of cluster id, group label and QC metrics (tsv)
%
%   See also LoadSpikeData, SpikeVisualizationApp.

    arguments
        phyDir (1,1) string {mustBeFolder}
        options.DataChannel (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        options.PreSamples (1,1) double {mustBeInteger, mustBePositive} = 40
        options.PostSamples (1,1) double {mustBeInteger, mustBePositive} = 40
        options.ExtractWaveforms (1,1) logical = true
    end

    phyDir = string(what(phyDir).path);   % resolve to an absolute path
    paramsFile = fullfile(phyDir, "params.py");
    if ~isfile(paramsFile)
        error("LoadSpikesPhy:noParams", ...
            "No params.py in %s. Point to a Phy export folder.", phyDir);
    end

    % --- recording parameters (parsed from params.py) ---
    params = readPhyParams(paramsFile);
    spikes.sourceType = "phy";
    spikes.phyDir = phyDir;
    spikes.datPath = fullfile(phyDir, params.datPath);
    spikes.samplingRate = params.sampleRate;
    spikes.gainToUV = params.gain;
    spikes.numChannels = params.numChannels;
    spikes.dataChannel = options.DataChannel;

    % --- per-spike arrays ---
    spikes.spikeTimes = double(readNPY(fullfile(phyDir, "spike_times.npy")));
    spikes.clusters = double(readNPY(fullfile(phyDir, "spike_clusters.npy")));
    spikes.templates = double(readNPY(fullfile(phyDir, "spike_templates.npy")));
    spikes.amplitudes = double(readNPY(fullfile(phyDir, "amplitudes.npy")));

    % Sort by time so downstream ISI/CCG math is well defined.
    [spikes.spikeTimes, order] = sort(spikes.spikeTimes);
    spikes.clusters = spikes.clusters(order);
    spikes.templates = spikes.templates(order);
    spikes.amplitudes = spikes.amplitudes(order);

    % --- sorter PC features (kept for the FeatureView / split tool) ---
    pcFile = fullfile(phyDir, "pc_features.npy");
    if isfile(pcFile)
        spikes.pcFeatures = single(readNPY(pcFile));
        spikes.pcFeatures = spikes.pcFeatures(order, :, :);
        spikes.pcFeatureInd = double(readNPY(fullfile(phyDir, "pc_feature_ind.npy")));
    else
        spikes.pcFeatures = single.empty(numel(spikes.spikeTimes), 0, 0);
        spikes.pcFeatureInd = [];
    end

    % --- raw data + waveform extraction ---
    numSamples = NaN;
    if isfile(spikes.datPath)
        totalInt16 = dir(spikes.datPath).bytes / 2;   % int16 = 2 bytes
        % Some exports declare the wrong n_channels_dat in params.py (e.g. a
        % 2-channel ch0-real/ch1-dup file labelled 1 channel). Reading it with
        % the wrong count interleaves the channels, so waveforms turn to noise
        % and spike markers land at the wrong time. Infer the true count from
        % where the spikes fall and how sharply they align.
        spikes.numChannels = inferChannelCount(spikes.datPath, totalInt16, ...
            spikes.numChannels, spikes.spikeTimes, spikes.clusters);
        numSamples = floor(totalInt16 / spikes.numChannels);
    end
    spikes.numSamples = numSamples;
    spikes.waveformWindow = [options.PreSamples options.PostSamples];

    if options.ExtractWaveforms && isfile(spikes.datPath) && ~isnan(numSamples)
        spikes.waveforms = extractWaveforms(spikes.datPath, numSamples, ...
            spikes.numChannels, options.DataChannel, spikes.spikeTimes, ...
            options.PreSamples, options.PostSamples);
    else
        spikes.waveforms = single.empty(numel(spikes.spikeTimes), 0);
    end

    % --- cluster metadata (group + Bombcell/QC metrics from tsv files) ---
    spikes.clusterTable = readClusterTables(phyDir, unique(spikes.clusters));
end

function params = readPhyParams(paramsFile)
%READPHYPARAMS Parse the handful of assignments used from a Phy params.py.
    txt = string(fileread(paramsFile));
    params.datPath = readValue(txt, "dat_path", "data.dat", true);
    params.numChannels = str2double(readValue(txt, "n_channels_dat", "1", false));
    params.sampleRate = str2double(readValue(txt, "sample_rate", "30000", false));
    params.gain = str2double(readValue(txt, "gain", "1", false));
    if isnan(params.gain)
        params.gain = 1;
    end
end

function value = readValue(txt, key, default, isString)
%READVALUE Return the right-hand side of `key = value` in a params.py text.
    pattern = key + "\s*=\s*([^\r\n]+)";   % value is the rest of that line only
    tok = regexp(txt, pattern, "tokens", "once");
    if isempty(tok)
        value = default;
        return
    end
    value = strtrim(tok(1));
    if isString
        value = erase(value, ["'" '"']);   % strip surrounding quotes
    end
end

function waveforms = extractWaveforms(datPath, numSamples, numChannels, ...
        dataChannel, spikeTimes, preSamples, postSamples)
%EXTRACTWAVEFORMS Pull one waveform per spike from a channel of a flat int16 file.
%   Spikes whose window would fall outside the recording are zero-padded.
    map = memmapfile(datPath, Format={'int16', [numChannels numSamples], 'raw'});
    channelData = single(map.Data.raw(dataChannel + 1, :));   % 1-by-numSamples

    windowLength = preSamples + postSamples + 1;
    numSpikes = numel(spikeTimes);
    waveforms = zeros(numSpikes, windowLength, "single");

    % spike_times.npy holds 0-BASED sample indices (Phy/KiloSort), so
    % sample t lives at MATLAB index t+1. Without the +1 every waveform
    % is cut one sample early, which shifts the template off the trough
    % and costs sensitivity wherever waveforms are matched.
    starts = round(spikeTimes) - preSamples + 1;
    for k = 1:numSpikes
        idx = starts(k):(starts(k) + windowLength - 1);
        valid = idx >= 1 & idx <= numSamples;
        waveforms(k, valid) = channelData(idx(valid));
    end
end

function nCh = inferChannelCount(datPath, totalInt16, nChParam, spikeTimes, clusters)
%INFERCHANNELCOUNT Best channel count for a flat int16 file, cross-checked
%   against the spikes. Trusts params unless another count both fits the spike
%   range and aligns the waveforms far more sharply (params can be wrong).
    nCh = nChParam;
    maxSt = max(spikeTimes);
    if isempty(maxSt) || maxSt <= 0
        return
    end
    % Candidate counts whose per-channel length still contains every spike.
    cand = unique([nChParam, 1, 2, round(totalInt16 / maxSt)]);
    cand = cand(cand >= 1 & cand <= 8 & floor(totalInt16 ./ cand) >= maxSt);
    if numel(cand) <= 1
        if ~isempty(cand)
            nCh = cand(1);
        end
        return
    end
    % Pick the count under which the dominant cluster's troughs line up best
    % (a wrong count interleaves channels -> troughs scatter at random phases).
    c0 = mode(clusters);
    s0 = spikeTimes(clusters == c0);
    bestScore = -inf;
    for k = 1:numel(cand)
        ns = floor(totalInt16 / cand(k));
        score = alignmentScore(datPath, cand(k), ns, s0);
        if score > bestScore
            bestScore = score;
            nCh = cand(k);
        end
    end
end

function score = alignmentScore(datPath, nCh, ns, times)
%ALIGNMENTSCORE Fraction of waveform troughs within +/-2 samples of their
%   median position (high for a correctly-read channel, low for interleaved).
    score = -inf;
    times = times(times > 50 & times < ns - 51);
    if numel(times) < 20
        return
    end
    times = times(round(linspace(1, numel(times), min(200, numel(times)))));
    map = memmapfile(datPath, Format={'int16', [nCh ns], 'raw'});
    troughs = zeros(numel(times), 1);
    for i = 1:numel(times)
        w = single(map.Data.raw(1, times(i) - 40 + 1:times(i) + 40 + 1));
        [~, troughs(i)] = min(w);
    end
    score = mean(abs(troughs - median(troughs)) <= 2);
end

function clusterTable = readClusterTables(phyDir, clusterIds)
%READCLUSTERTABLES Merge every cluster_*.tsv into one table keyed by cluster id.
    clusterTable = table(clusterIds(:), VariableNames="cluster_id");
    tsvFiles = dir(fullfile(phyDir, "cluster_*.tsv"));
    for k = 1:numel(tsvFiles)
        tsvPath = fullfile(tsvFiles(k).folder, tsvFiles(k).name);
        try
            t = readtable(tsvPath, FileType="text", Delimiter="\t", ...
                VariableNamingRule="preserve");
        catch
            continue   % skip malformed / empty tsv rather than failing the load
        end
        if ~ismember("cluster_id", string(t.Properties.VariableNames)) || width(t) < 2
            continue
        end
        clusterTable = outerjoin(clusterTable, t, Keys="cluster_id", ...
            MergeKeys=true, Type="left");
    end
end
