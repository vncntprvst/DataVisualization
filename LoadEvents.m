function ev = LoadEvents(source)
%LOADEVENTS Load event alignment times (recording seconds) for PETHs.
%   ev = LOADEVENTS(folder) looks for events.mat (then events.tsv / events.csv)
%   in the folder. ev = LOADEVENTS(file) loads that file directly.
%
%   Returns a struct:
%       ev.events  - struct: field per event name -> row vector of alignment
%                    times in recording seconds (same clock as spike_times/fs).
%       ev.cond    - struct: field per event name -> string array of condition
%                    labels (one per alignment time), where available.
%       ev.meta    - struct of metadata (as saved), or struct() if none.
%       ev.names   - string array of event names.
%       ev.source  - resolved file path.
%
%   Supported files:
%     * .mat with a struct `events` (field -> times, seconds); optional
%       `events_cond` (field -> conditions) and `meta`.
%     * .tsv / .csv in long form with columns `event`, `time_s`, and optional
%       `condition` (one row per event occurrence).
%
%   See also PETHTool, SpikeVisualizationApp.

    arguments
        source (1,1) string
    end

    file = resolveFile(source);
    if file == ""
        error("LoadEvents:notFound", ...
            "No events.mat/.tsv/.csv found for %s.", source);
    end

    [~, ~, ext] = fileparts(file);
    switch lower(ext)
        case ".mat"
            ev = fromMat(file);
        case {".tsv", ".csv"}
            ev = fromTable(file, ext);
        otherwise
            error("LoadEvents:badType", "Unsupported events file: %s", file);
    end
    ev.source = file;
    ev.names = string(fieldnames(ev.events));
end

function file = resolveFile(source)
    file = "";
    if isfolder(source)
        for name = ["events.mat", "events.tsv", "events.csv"]
            candidate = fullfile(source, name);
            if isfile(candidate)
                file = string(candidate);
                return
            end
        end
    elseif isfile(source)
        file = source;
    end
end

function ev = fromMat(file)
    s = load(file);
    if ~isfield(s, "events") || ~isstruct(s.events)
        error("LoadEvents:noEvents", ...
            "%s has no `events` struct (field -> times in seconds).", file);
    end
    ev.events = structfun(@(x) double(x(:)'), s.events, UniformOutput=false);
    ev.cond = struct();
    if isfield(s, "events_cond") && isstruct(s.events_cond)
        ev.cond = structfun(@(x) string(x(:)'), s.events_cond, ...
            UniformOutput=false);
    end
    ev.meta = struct();
    if isfield(s, "meta") && isstruct(s.meta)
        ev.meta = s.meta;
    end
end

function ev = fromTable(file, ext)
    delim = "\t";
    if lower(ext) == ".csv"
        delim = ",";
    end
    t = readtable(file, FileType="text", Delimiter=delim, ...
        TextType="string", VariableNamingRule="preserve");
    vars = string(t.Properties.VariableNames);
    if ~all(ismember(["event", "time_s"], vars))
        error("LoadEvents:badColumns", ...
            "%s needs columns 'event' and 'time_s'.", file);
    end
    ev.events = struct();
    ev.cond = struct();
    ev.meta = struct();
    hasCond = ismember("condition", vars);
    names = unique(t.event, "stable");
    for name = names'
        rows = t.event == name;
        fieldName = matlab.lang.makeValidName(name);
        ev.events.(fieldName) = double(t.time_s(rows))';
        if hasCond
            ev.cond.(fieldName) = string(t.condition(rows))';
        end
    end
end
