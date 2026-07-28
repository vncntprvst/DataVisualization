function ds = LoadDataset(source)
%LOADDATASET Load a dataset table (one row per unit) for the DatasetBrowser.
%   ds = LOADDATASET(file) reads a tab- or comma-separated table describing a
%   set of sorted units / recordings (e.g. D:\CB paper\0-data\dataset.tsv) and
%   returns a struct the browser can display and filter.
%
%   Returns a struct:
%       ds.table            the raw table (one row per unit).
%       ds.vars             string array of column names.
%       ds.source           resolved file path.
%       ds.hasRecordingTag  true if a `recording_tag` column is present (needed
%                           to locate the Phy sorting folder to open).
%
%   Key text columns (recording_tag, unit, region, task, bombcell_label) are
%   forced to text so integer-looking tags/labels are not coerced to numbers.
%
%   See also DatasetBrowser, LoadEvents, SpikeVisualizationApp.

    arguments
        source (1,1) string {mustBeFile}
    end

    [~, ~, ext] = fileparts(source);
    opts = detectImportOptions(source, FileType="text", ...
        VariableNamingRule="preserve");
    if lower(ext) == ".tsv"
        opts.Delimiter = "\t";
    end

    % Keep identifier-like columns as text (tags such as S113L4A5_12540 and
    % labels such as neg1+pos4+pos6 must not become numbers / NaN).
    textCols = ["recording_tag", "unit", "region", "task", "bombcell_label", ...
        "sorter", "is_somatic", "purkinje_flag", "in_original", "orig_basis"];
    present = intersect(textCols, string(opts.VariableNames));
    if ~isempty(present)
        opts = setvartype(opts, cellstr(present), "string");
    end

    t = readtable(source, opts);
    ds.table = t;
    ds.vars = string(t.Properties.VariableNames);
    ds.source = string(source);
    ds.hasRecordingTag = ismember("recording_tag", ds.vars);
end
