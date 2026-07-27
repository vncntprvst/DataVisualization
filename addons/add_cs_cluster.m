function csClusterId = add_cs_cluster(phyDir, options)
% AddonName: Add CS cluster
%ADD_CS_CLUSTER Inject the P-sort complex-spike train into a Phy folder.
%   add_cs_cluster(phyDir) finds the session's P-sort Purkinje output
%   (purkinje/purkinje_*.npz, searched in phyDir and its parents), and appends
%   the complex-spike train (cs_train) to the Phy sorting as its own new
%   cluster. Any Phy-aware tool - Phy itself or SpikeVisualizationApp - then
%   shows the complex spikes as a normal, separate cluster. The two trains are
%   NOT pooled: the CS spikes get a distinct cluster id, never merged into a
%   simple-spike cluster.
%
%   csClusterId = add_cs_cluster(phyDir, Name=Value) accepts:
%       PurkinjeNpz - path to the .npz to use (default: auto-find).
%       Force       - re-inject even if a previous injection is detected.
%
%   The per-spike Phy arrays (spike_times, spike_clusters, spike_templates,
%   amplitudes, pc_features) are extended, re-sorted by time and rewritten; the
%   originals are backed up once as *.preCS.npy. cluster_group.tsv gets a "good"
%   row for the new cluster, and a .cs_injected.json marker records the id (used
%   to stay idempotent). To undo, restore the *.preCS.npy files.
%
%   Requires npy-matlab (readNPY/writeNPY) on the path.
%
%   See also LoadSpikesPhy, SpikeVisualizationApp.

    arguments
        phyDir (1,1) string {mustBeFolder}
        options.PurkinjeNpz (1,1) string = ""
        options.Force (1,1) logical = false
    end

    phyDir = string(what(phyDir).path);
    marker = fullfile(phyDir, ".cs_injected.json");
    if isfile(marker) && ~options.Force
        info = jsondecode(fileread(marker));
        csClusterId = info.csClusterId;
        fprintf("CS already injected as cluster %d (use Force=true to redo).\n", ...
            csClusterId);
        return
    end

    npz = options.PurkinjeNpz;
    if npz == ""
        npz = findPurkinjeNpz(phyDir);
    end
    if npz == "" || ~isfile(npz)
        error("add_cs_cluster:noPurkinje", ...
            "No purkinje_*.npz found for %s. Pass PurkinjeNpz explicitly.", phyDir);
    end

    % --- read the CS train (samples) ---
    cs = readPurkinjeCS(npz);
    fs = readPhySampleRate(fullfile(phyDir, "params.py"));
    numSamples = countSamples(phyDir);
    csTimes = cs.times(:);
    if ~isnan(numSamples) && max(csTimes) <= (numSamples / fs) * 2
        csTimes = csTimes * fs;         % stored in seconds after all
    end
    csTimes = round(csTimes);
    nCS = numel(csTimes);

    % --- read existing per-spike arrays ---
    spikeTimes = int64(readNPY(fullfile(phyDir, "spike_times.npy")));
    clusters = int32(readNPY(fullfile(phyDir, "spike_clusters.npy")));
    templates = int32(readNPY(fullfile(phyDir, "spike_templates.npy")));
    amplitudes = double(readNPY(fullfile(phyDir, "amplitudes.npy")));

    csClusterId = double(max(clusters)) + 1;
    csTemplate = int32(0);              % reuse template 0 to stay Phy-valid
    csAmplitude = 1;
    if ~isempty(cs.waveform)
        csAmplitude = max(cs.waveform) - min(cs.waveform);
    end

    spikeTimes = [spikeTimes; int64(csTimes)];
    clusters = [clusters; repmat(int32(csClusterId), nCS, 1)];
    templates = [templates; repmat(csTemplate, nCS, 1)];
    amplitudes = [amplitudes; repmat(csAmplitude, nCS, 1)];

    pcFile = fullfile(phyDir, "pc_features.npy");
    hasPc = isfile(pcFile);
    if hasPc
        pcFeatures = single(readNPY(pcFile));
        pad = zeros([nCS, size(pcFeatures, 2), size(pcFeatures, 3)], "single");
        pcFeatures = cat(1, pcFeatures, pad);
    end

    % --- keep everything time-sorted (Phy expects ascending spike_times) ---
    [spikeTimes, order] = sort(spikeTimes);
    clusters = clusters(order);
    templates = templates(order);
    amplitudes = amplitudes(order);
    if hasPc
        pcFeatures = pcFeatures(order, :, :);
    end

    % --- write atomically: all to temp files first, then move into place ---
    % This keeps the folder consistent even if a write fails partway (e.g. a
    % cloud-sync lock or a read-only file), rolling back any files moved.
    targets = ["spike_times.npy", "spike_clusters.npy", ...
        "spike_templates.npy", "amplitudes.npy"];
    payload = {spikeTimes, clusters, templates, amplitudes};
    if hasPc
        targets(end + 1) = "pc_features.npy";
        payload{end + 1} = pcFeatures;
    end
    tmpFiles = fullfile(phyDir, targets + ".cstmp");

    try
        for i = 1:numel(targets)
            writeNPY(payload{i}, tmpFiles(i));
        end
    catch err
        deleteIfExist(tmpFiles);
        error("add_cs_cluster:writeTemp", ...
            "Could not write temporary files in %s: %s", phyDir, err.message);
    end

    moved = strings(1, 0);
    try
        for i = 1:numel(targets)
            backupOnce(phyDir, targets(i));
            moveWithRetry(tmpFiles(i), fullfile(phyDir, targets(i)));
            moved(end + 1) = targets(i); %#ok<AGROW>
        end
    catch err
        for f = moved   % roll back the files already replaced
            backup = fullfile(phyDir, replace(f, ".npy", ".preCS.npy"));
            if isfile(backup)
                ensureWritable(fullfile(phyDir, f));
                copyfile(backup, fullfile(phyDir, f));
            end
        end
        deleteIfExist(tmpFiles);
        error("add_cs_cluster:move", ...
            "A sorting file is locked by another program, so nothing was " + ...
            "changed (folder left intact). Close Phy and any file-sync/antivirus " + ...
            "using this folder, then retry. Details: %s", err.message);
    end

    addClusterGroupRow(phyDir, csClusterId, "good");

    info = struct("csClusterId", csClusterId, "nCS", nCS, ...
        "source", npz, "label", "CS");
    try
        ensureWritable(marker);
        writelines(jsonencode(info), marker);
    catch
        % the marker is only used for idempotency; injection already succeeded
    end

    fprintf("Injected %d complex spikes as cluster %d in %s\n", ...
        nCS, csClusterId, phyDir);
end

% =========================================================================
function npz = findPurkinjeNpz(phyDir)
    npz = "";
    d = char(phyDir);
    for up = 1:4
        hits = dir(fullfile(d, "purkinje", "purkinje_*.npz"));
        if ~isempty(hits)
            npz = string(fullfile(hits(1).folder, hits(1).name));
            return
        end
        parent = fileparts(d);
        if strcmp(parent, d)
            return
        end
        d = parent;
    end
end

function cs = readPurkinjeCS(npzFile)
    tmp = tempname;
    mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, "s"));
    unzip(npzFile, tmp);
    cs.times = double(readNPY(fullfile(tmp, "cs_train.npy")));
    waveformFile = fullfile(tmp, "cs_waveform.npy");
    if isfile(waveformFile)
        cs.waveform = double(readNPY(waveformFile));
    else
        cs.waveform = [];
    end
    jsonFile = replace(npzFile, ".npz", ".json");
    if isfile(jsonFile)
        cs.meta = jsondecode(fileread(jsonFile));
    else
        cs.meta = struct();
    end
end

function fs = readPhySampleRate(paramsFile)
    fs = 30000;
    if ~isfile(paramsFile)
        return
    end
    txt = string(fileread(paramsFile));
    tok = regexp(txt, "sample_rate\s*=\s*([^\r\n]+)", "tokens", "once");
    if ~isempty(tok)
        fs = str2double(tok(1));
    end
end

function numSamples = countSamples(phyDir)
    numSamples = NaN;
    datFile = fullfile(phyDir, "data.dat");
    paramsFile = fullfile(phyDir, "params.py");
    if ~isfile(datFile) || ~isfile(paramsFile)
        return
    end
    txt = string(fileread(paramsFile));
    tok = regexp(txt, "n_channels_dat\s*=\s*([^\r\n]+)", "tokens", "once");
    numChannels = 1;
    if ~isempty(tok)
        numChannels = str2double(tok(1));
    end
    numSamples = dir(datFile).bytes / 2 / numChannels;   % int16 = 2 bytes
end

function backupOnce(phyDir, fileName)
    [~, base, ext] = fileparts(fileName);
    src = fullfile(phyDir, fileName);
    backup = fullfile(phyDir, base + ".preCS" + ext);
    if isfile(src) && ~isfile(backup)
        copyfile(src, backup);
    end
end

function moveWithRetry(src, dst)
    %MOVEWITHRETRY Move src onto dst, retrying through transient file locks
    %   (e.g. file-sync or antivirus briefly holding the target open).
    ensureWritable(dst);
    lastErr = [];
    for attempt = 1:6
        try
            movefile(src, dst);
            return
        catch err
            lastErr = err;
            pause(0.4 * attempt);
            ensureWritable(dst);
        end
    end
    rethrow(lastErr);
end

function ensureWritable(file)
    %ENSUREWRITABLE Clear a read-only attribute so the file can be overwritten.
    if isfile(file)
        try
            fileattrib(file, "+w");
        catch
            % best effort; the caller handles a genuine write failure
        end
    end
end

function deleteIfExist(files)
    for f = files
        if isfile(f)
            try
                delete(f);
            catch
                % leftover temp file is harmless
            end
        end
    end
end

function addClusterGroupRow(phyDir, clusterId, group)
    file = fullfile(phyDir, "cluster_group.tsv");
    if isfile(file)
        t = readtable(file, FileType="text", Delimiter="\t", ...
            VariableNamingRule="preserve");
    else
        t = table('Size', [0 2], VariableTypes=["double", "string"], ...
            VariableNames=["cluster_id", "group"]);
    end
    if ~any(t.cluster_id == clusterId)
        t = [t; {clusterId, string(group)}];
    end
    ensureWritable(file);
    writetable(t, file, FileType="text", Delimiter="\t");
end
