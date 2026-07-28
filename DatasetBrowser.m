classdef DatasetBrowser < handle
    %DATASETBROWSER Browse a dataset table and open its sortings/units.
    %   DatasetBrowser opens the browser as a standalone window; the first
    %   time you open an entry it launches a SpikeVisualizationApp for it (and
    %   reuses that window afterwards).
    %
    %   DatasetBrowser(app) attaches to an existing SpikeVisualizationApp
    %   (this is how the app's Tools menu opens it).
    %
    %   The browser reads a dataset table (one row per unit, e.g.
    %   D:\CB paper\0-data\dataset.tsv) and lets you browse it by Units or
    %   Recordings, filter and sort, then open an entry in the app.
    %
    %   Opening a Recording loads its Phy sorting (<sortings root>\<recording_tag>
    %   \phy_postmerge). Opening a Unit loads that sorting, selects the unit's
    %   cluster (resolved via final_labels.tsv) and opens the CurationTool on it.
    %   Tags with no Phy folder fall back to their online-sorted (REX) session,
    %   loaded spikes-only.
    %
    %   See also LoadDataset, LoadSpikesOnline, SpikeVisualizationApp.

    properties (Access = private)
        App                        % SpikeVisualizationApp handle
        Dataset struct = struct()  % from LoadDataset
        RawTable table             % one row per unit
        ViewTable table            % table currently shown (filtered + granularity)
        RowKind string = "units"   % "units" | "recordings"
        SortingsRoot string = ""   % <root>\<recording_tag>\phy_postmerge
        LoadedTag string = ""      % recording_tag currently loaded (for highlight)

        UIFigure       matlab.ui.Figure
        DatasetLabel   matlab.ui.control.Label
        RootLabel      matlab.ui.control.Label
        GranDropdown   matlab.ui.control.DropDown
        FilterDrops    struct = struct()   % column -> dropdown handle
        SearchField    matlab.ui.control.EditField
        Table          matlab.ui.control.Table
        DetailArea     matlab.ui.control.TextArea
        OpenButton     matlab.ui.control.Button
        InfoLabel      matlab.ui.control.Label
    end

    properties (Constant, Access = private)
        PrefGroup = "SpikeVisualizationApp"
        DefaultRoot = "D:\CB paper\1-analyses\spikesort_cluster\curation_all"
        DefaultDataset = "D:\CB paper\0-data\dataset.tsv"
        ProcessedRoot = "D:\CB paper\0-data\SLData\processed"   % REX online sessions
        % Filter bar: {display label, column name}.
        FilterCols = ["region", "task", "bombcell_label", "purkinje_flag"]
        % Preferred column order for the Units view (intersected with what's there).
        UnitsPreset = ["recording_tag", "unit", "region", "task", ...
            "bombcell_label", "single_unit", "snr", "isi_viol_2p0_pct", ...
            "presence_ratio", "n_spikes", "purkinje_flag", "cs_rate_hz"]
    end

    methods
        function tool = DatasetBrowser(app)
            arguments
                app = SpikeVisualizationApp.empty   % [] -> standalone; created on first open
            end
            if ~isempty(app) && isa(app, "SpikeVisualizationApp") && isvalid(app)
                tool.App = app;
            else
                tool.App = SpikeVisualizationApp.empty;
            end
            tool.buildUI();
            tool.resolveRoot();          % confirm sortings root on first use
            tool.syncLoadedTagFromApp();
            startPath = string(getpref(tool.PrefGroup, "datasetPath", ...
                tool.DefaultDataset));
            tool.tryLoadDataset(startPath);
        end

        function delete(tool)
            if ~isempty(tool.UIFigure) && isvalid(tool.UIFigure)
                delete(tool.UIFigure);
            end
        end

        function f = figureHandle(tool)
            f = tool.UIFigure;
        end
    end

    % ------------------------------------------------------------------- UI
    methods (Access = private)
        function buildUI(tool)
            tool.UIFigure = uifigure(Name="Dataset browser", ...
                Position=[140 130 1120 720]);
            outer = uigridlayout(tool.UIFigure, [3 1]);
            outer.RowHeight = {"fit", "1x", 22};
            outer.Padding = [6 6 6 6];

            % --- top controls (config bar + filter bar) ---
            controls = uigridlayout(outer, [2 1]);
            controls.Layout.Row = 1;
            controls.RowHeight = {30, 30};
            controls.Padding = [0 0 0 0];
            controls.RowSpacing = 4;

            barA = uigridlayout(controls, [1 6]);
            barA.Layout.Row = 1;
            barA.ColumnWidth = {120, "1x", 120, "1x", 80, 130};
            barA.Padding = [0 0 0 0];
            uibutton(barA, Text="Load dataset...", ...
                ButtonPushedFcn=@(~, ~) tool.loadDatasetDialog());
            tool.DatasetLabel = uilabel(barA, Text="(no dataset)", ...
                FontColor=[0.3 0.3 0.3]);
            uibutton(barA, Text="Sortings root...", ...
                ButtonPushedFcn=@(~, ~) tool.chooseRootDialog());
            tool.RootLabel = uilabel(barA, Text="", FontColor=[0.3 0.3 0.3]);
            uilabel(barA, Text="Granularity", HorizontalAlignment="right");
            tool.GranDropdown = uidropdown(barA, ...
                Items=["Units", "Recordings"], ItemsData=["units", "recordings"], ...
                Value="units", ValueChangedFcn=@(s, ~) tool.setGranularity(s.Value));

            barB = uigridlayout(controls, [1 11]);
            barB.Layout.Row = 2;
            barB.ColumnWidth = {"fit", 110, "fit", 110, "fit", 110, "fit", 110, ...
                "fit", "1x", 60};
            barB.Padding = [0 0 0 0];
            for k = 1:numel(tool.FilterCols)
                col = tool.FilterCols(k);
                uilabel(barB, Text=filterLabel(col), HorizontalAlignment="right");
                d = uidropdown(barB, Items={'(all)'}, ...
                    ValueChangedFcn=@(~, ~) tool.applyFilters());
                tool.FilterDrops.(col) = d;
            end
            uilabel(barB, Text="Search", HorizontalAlignment="right");
            tool.SearchField = uieditfield(barB, "text", ...
                ValueChangedFcn=@(~, ~) tool.applyFilters());
            uibutton(barB, Text="Clear", ...
                ButtonPushedFcn=@(~, ~) tool.clearFilters());

            % --- table + detail pane ---
            mid = uigridlayout(outer, [1 2]);
            mid.Layout.Row = 2;
            mid.ColumnWidth = {"1x", 300};
            mid.Padding = [0 0 0 0];

            tool.Table = uitable(mid, ColumnSortable=true, ...
                SelectionType="row", Multiselect="off", RowName={}, ...
                SelectionChangedFcn=@(s, e) tool.onRowSelected(e), ...
                DoubleClickedFcn=@(s, e) tool.onOpen());

            right = uigridlayout(mid, [3 1]);
            right.RowHeight = {"fit", "1x", 30};
            right.Padding = [0 0 0 0];
            uilabel(right, Text="Selected entry", FontWeight="bold");
            tool.DetailArea = uitextarea(right, Editable="off", Value="");
            tool.OpenButton = uibutton(right, Text="Open selected", ...
                BackgroundColor=[0.85 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.onOpen());

            tool.InfoLabel = uilabel(outer, Text="");
            tool.InfoLabel.Layout.Row = 3;
        end
    end

    % ------------------------------------------------------- root / dataset
    methods (Access = private)
        function resolveRoot(tool)
            %RESOLVEROOT Use the stored sortings root, else confirm one with the
            %   user on first use and remember the choice.
            root = string(getpref(tool.PrefGroup, "sortingsRoot", ""));
            if root ~= "" && isfolder(root)
                tool.setRoot(root, false);
                return
            end
            cand = tool.DefaultRoot;
            if ~isfolder(cand)
                cand = string(pwd);
            end
            sel = uigetdir(char(cand), ...
                "Confirm the sortings root (contains <recording_tag>\phy_postmerge)");
            figure(tool.UIFigure);
            if isequal(sel, 0)
                tool.setRoot(cand, false);   % use default without persisting
            else
                tool.setRoot(string(sel), true);
            end
        end

        function chooseRootDialog(tool)
            start = tool.SortingsRoot;
            if start == "" || ~isfolder(start)
                start = tool.DefaultRoot;
            end
            sel = uigetdir(char(start), "Select the sortings root");
            figure(tool.UIFigure);
            if isequal(sel, 0)
                return
            end
            tool.setRoot(string(sel), true);
            tool.renderTable();   % refresh missing-folder styling
        end

        function setRoot(tool, root, persist)
            tool.SortingsRoot = root;
            tool.RootLabel.Text = compactPath(root);
            tool.RootLabel.Tooltip = root;
            if persist
                setpref(tool.PrefGroup, "sortingsRoot", char(root));
            end
        end

        function tryLoadDataset(tool, source)
            if source == "" || ~isfile(source)
                tool.setInfo("Load a dataset with 'Load dataset...' " + ...
                    "(default not found: " + source + ").");
                return
            end
            try
                tool.setDataset(LoadDataset(source));
                setpref(tool.PrefGroup, "datasetPath", char(source));
            catch err
                tool.setInfo("Dataset load failed: " + err.message);
            end
        end

        function loadDatasetDialog(tool)
            [f, d] = uigetfile({'*.tsv;*.csv', "Dataset files"}, ...
                "Select a dataset table");
            figure(tool.UIFigure);
            if isequal(f, 0)
                return
            end
            tool.tryLoadDataset(string(fullfile(d, f)));
        end

        function setDataset(tool, ds)
            tool.Dataset = ds;
            tool.RawTable = ds.table;
            tool.DatasetLabel.Text = compactPath(ds.source);
            tool.DatasetLabel.Tooltip = ds.source;
            if ~ds.hasRecordingTag
                tool.OpenButton.Enable = "off";
                tool.GranDropdown.Value = "units";
                tool.GranDropdown.Enable = "off";
                tool.RowKind = "units";
            else
                tool.OpenButton.Enable = "on";
                tool.GranDropdown.Enable = "on";
            end
            tool.populateFilters();
            tool.applyFilters();
        end

        function populateFilters(tool)
            for k = 1:numel(tool.FilterCols)
                col = tool.FilterCols(k);
                d = tool.FilterDrops.(col);
                if ismember(col, tool.Dataset.vars)
                    vals = unique(string(tool.RawTable.(col)), "stable");
                    vals = sort(vals(~ismissing(vals) & vals ~= ""));
                    d.Items = ["(all)"; vals(:)]';
                    d.Value = "(all)";
                    d.Enable = "on";
                else
                    d.Items = {'(all)'};
                    d.Value = "(all)";
                    d.Enable = "off";
                end
            end
        end
    end

    % ----------------------------------------------------- filter / display
    methods (Access = private)
        function setGranularity(tool, kind)
            tool.RowKind = string(kind);
            tool.applyFilters();
        end

        function clearFilters(tool)
            for k = 1:numel(tool.FilterCols)
                d = tool.FilterDrops.(tool.FilterCols(k));
                if ~isempty(d.Items)
                    d.Value = "(all)";
                end
            end
            tool.SearchField.Value = "";
            tool.applyFilters();
        end

        function applyFilters(tool)
            if isempty(tool.RawTable)
                return
            end
            t = tool.RawTable;
            keep = true(height(t), 1);
            for k = 1:numel(tool.FilterCols)
                col = tool.FilterCols(k);
                if ~ismember(col, tool.Dataset.vars)
                    continue
                end
                val = string(tool.FilterDrops.(col).Value);
                if val ~= "(all)"
                    keep = keep & (string(t.(col)) == val);
                end
            end
            q = lower(strtrim(string(tool.SearchField.Value)));
            if q ~= ""
                hit = false(height(t), 1);
                for col = ["recording_tag", "unit"]
                    if ismember(col, tool.Dataset.vars)
                        hit = hit | contains(lower(string(t.(col))), q);
                    end
                end
                keep = keep & hit;
            end
            filtered = t(keep, :);

            if tool.RowKind == "recordings" && tool.Dataset.hasRecordingTag
                tool.ViewTable = tool.aggregateByRecording(filtered);
            else
                tool.ViewTable = filtered;
            end
            tool.renderTable();
        end

        function agg = aggregateByRecording(~, t)
            if isempty(t) || ~ismember("recording_tag", string(t.Properties.VariableNames))
                agg = table();
                return
            end
            tags = unique(string(t.recording_tag), "stable");
            n = numel(tags);
            nUnits = zeros(n, 1);
            region = strings(n, 1);
            task = strings(n, 1);
            for i = 1:n
                rows = string(t.recording_tag) == tags(i);
                nUnits(i) = sum(rows);
                region(i) = firstValue(t, "region", rows);
                task(i) = firstValue(t, "task", rows);
            end
            agg = table(tags, nUnits, region, task, ...
                VariableNames=["recording_tag", "n_units", "region", "task"]);
        end

        function vars = displayVars(tool)
            if tool.RowKind == "recordings"
                vars = string(tool.ViewTable.Properties.VariableNames);
            else
                vars = intersect(tool.UnitsPreset, tool.Dataset.vars, "stable");
                if isempty(vars)
                    vars = tool.Dataset.vars;
                end
            end
        end

        function renderTable(tool)
            if isempty(tool.ViewTable)
                tool.Table.Data = table();
                tool.setInfo("No rows match the current filters.");
                return
            end
            vars = tool.displayVars();
            tool.Table.Data = tool.ViewTable(:, vars);
            removeStyle(tool.Table);

            tags = string(tool.ViewTable.recording_tag);
            % Grey out rows whose sorting folder (and online fallback) is absent.
            missing = false(numel(tags), 1);
            for i = 1:numel(tags)
                missing(i) = ~tool.hasSorting(tags(i));
            end
            if any(missing)
                addStyle(tool.Table, uistyle(FontColor=[0.6 0.6 0.6]), ...
                    "row", find(missing));
            end
            % Highlight the recording currently loaded in the app.
            if tool.LoadedTag ~= ""
                hit = find(tags == tool.LoadedTag);
                if ~isempty(hit)
                    addStyle(tool.Table, uistyle(BackgroundColor=[0.83 0.94 0.83]), ...
                        "row", hit);
                end
            end
            tool.setInfo(sprintf("%d rows shown (%s).", height(tool.ViewTable), ...
                tool.RowKind));
        end

        function onRowSelected(tool, ~)
            row = tool.selectedRow();
            if isempty(row)
                tool.DetailArea.Value = "";
                return
            end
            vars = string(tool.ViewTable.Properties.VariableNames);
            lines = strings(numel(vars), 1);
            for i = 1:numel(vars)
                lines(i) = vars(i) + ": " + string(tool.ViewTable.(vars(i))(row));
            end
            tool.DetailArea.Value = cellstr(lines);
            tag = string(tool.ViewTable.recording_tag(row));
            tool.OpenButton.Enable = ...
                matlab.lang.OnOffSwitchState(tool.hasSorting(tag));
        end
    end

    % ------------------------------------------------------------- opening
    methods (Access = private)
        function onOpen(tool)
            row = tool.selectedRow();
            if isempty(row)
                tool.setInfo("Select a row first.");
                return
            end
            tag = string(tool.ViewTable.recording_tag(row));
            if tag == "" || ismissing(tag)
                tool.setInfo("This row has no recording_tag.");
                return
            end
            phyDir = tool.phyDirFor(tag);
            try
                if isfolder(phyDir) && isfile(fullfile(phyDir, "params.py"))
                    if tool.appValid()
                        tool.App.loadPhyFolder(phyDir);
                    else
                        tool.discardStaleApp();
                        tool.App = SpikeVisualizationApp(phyDir);   % launch the GUI
                    end
                    mode = "phy";
                else
                    mat = tool.onlineMatFor(tag);
                    if mat == ""
                        tool.setInfo("No sorting found for " + tag + ...
                            " (looked in " + phyDir + " and the online sessions).");
                        return
                    end
                    model = LoadSpikesOnline(mat);
                    if tool.appValid()
                        tool.App.setData(model);
                    else
                        tool.discardStaleApp();
                        tool.App = SpikeVisualizationApp(model);    % launch the GUI
                    end
                    mode = "online";
                end
            catch err
                tool.setInfo("Open failed: " + err.message);
                return
            end

            tool.App.refreshCompanions();

            if mode == "phy" && tool.RowKind == "units" ...
                    && ismember("unit", string(tool.ViewTable.Properties.VariableNames))
                label = string(tool.ViewTable.unit(row));
                cid = tool.resolveClusterId(tag, label);
                if ~isnan(cid)
                    tool.App.setSelectedClusterIds(cid);
                    tool.App.openCuration(cid);   % skips silently if no waveforms
                else
                    tool.setInfo("Loaded " + tag + "; unit '" + label + ...
                        "' not found in final_labels.tsv (showing default cluster).");
                end
            end

            tool.LoadedTag = tag;
            tool.renderTable();
            if mode == "online"
                tool.setInfo("Loaded online (spikes-only) session for " + tag + ".");
            end
        end

        function tf = hasSorting(tool, tag)
            tf = isfolder(tool.phyDirFor(tag)) || tool.onlineMatFor(tag) ~= "";
        end

        function phyDir = phyDirFor(tool, tag)
            phyDir = fullfile(tool.SortingsRoot, tag, "phy_postmerge");
        end

        function mat = onlineMatFor(tool, tag)
            %ONLINEMATFOR Locate the REX online session file for a tag, if any.
            mat = "";
            if ~isfolder(tool.ProcessedRoot)
                return
            end
            hits = dir(fullfile(tool.ProcessedRoot, "*", tag + "_REX.mat"));
            if ~isempty(hits)
                mat = string(fullfile(hits(1).folder, hits(1).name));
            end
        end

        function cid = resolveClusterId(tool, tag, label)
            %RESOLVECLUSTERID Map a unit's merged_label to its integer cluster id
            %   via <root>\<tag>\final_labels.tsv (fallback: root\ALL_UNITS.tsv).
            cid = NaN;
            flPath = fullfile(tool.SortingsRoot, tag, "final_labels.tsv");
            if isfile(flPath)
                cid = lookupId(flPath, "final_id", "merged_label", label);
            end
            if isnan(cid)
                allPath = fullfile(tool.SortingsRoot, "ALL_UNITS.tsv");
                if isfile(allPath)
                    cid = lookupId(allPath, "final_unit_id", "merged_label", ...
                        label, "recording_tag", tag);
                end
            end
        end

        function row = selectedRow(tool)
            row = tool.Table.Selection;
            if isempty(row) || isempty(tool.ViewTable)
                row = [];
            else
                row = row(1);
            end
        end

        function tf = appValid(tool)
            % The app object can outlive its window (closing the figure does not
            % delete the handle), so also require a live figure — otherwise
            % loading into it plots onto deleted axes.
            tf = ~isempty(tool.App) && isvalid(tool.App) ...
                && isvalid(tool.App.figureHandle());
        end

        function discardStaleApp(tool)
            % Drop a defunct app (window closed) so a fresh one is launched.
            if ~isempty(tool.App) && isvalid(tool.App)
                delete(tool.App);
            end
            tool.App = SpikeVisualizationApp.empty;
        end

        function syncLoadedTagFromApp(tool)
            if ~tool.appValid()
                return
            end
            s = tool.App.Spikes;
            if isfield(s, "phyDir") && ~isempty(s.phyDir)
                parts = split(string(s.phyDir), filesep);
                if numel(parts) >= 2
                    tool.LoadedTag = parts(end - 1);   % <tag>\phy_postmerge
                end
            end
        end

        function setInfo(tool, msg)
            if ~isempty(tool.InfoLabel) && isvalid(tool.InfoLabel)
                tool.InfoLabel.Text = string(msg);
            end
        end
    end
end

% =========================================================================
% Local helper functions
% =========================================================================
function s = filterLabel(col)
    map = struct(region="Region", task="Task", ...
        bombcell_label="Bombcell", purkinje_flag="Purkinje");
    if isfield(map, col)
        s = map.(col);
    else
        s = string(col);
    end
end

function v = firstValue(t, col, rows)
    v = "";
    if ismember(col, string(t.Properties.VariableNames))
        vals = string(t.(col)(rows));
        if ~isempty(vals)
            v = vals(1);
        end
    end
end

function p = compactPath(fullPath)
    %COMPACTPATH Show the last two path components for a label.
    parts = split(string(fullPath), filesep);
    if numel(parts) >= 2
        p = "..." + filesep + join(parts(end - 1:end), filesep);
    else
        p = string(fullPath);
    end
end

function id = lookupId(tsvPath, idCol, labelCol, label, tagCol, tag)
    %LOOKUPID First idCol value whose labelCol matches label (and tagCol==tag).
    id = NaN;
    try
        t = readtable(tsvPath, FileType="text", Delimiter="\t", ...
            VariableNamingRule="preserve", TextType="string");
    catch
        return
    end
    vars = string(t.Properties.VariableNames);
    if ~all(ismember([idCol, labelCol], vars))
        return
    end
    match = string(t.(labelCol)) == label;
    if nargin >= 6 && ismember(tagCol, vars)
        match = match & string(t.(tagCol)) == tag;
    end
    hit = find(match, 1);
    if ~isempty(hit)
        id = double(t.(idCol)(hit));
    end
end
