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
        LoadedUnit string = ""     % unit label currently loaded (units view)
        LoadedLevel string = ""    % "unit" | "recording" of the loaded entry

        Locked logical = false     % when true, review state is read-only
        StateMap                   % containers.Map: key -> struct(status,note,modified)
        NoteField      matlab.ui.control.EditField
        LockBox        matlab.ui.control.CheckBox
        StatusButtons  matlab.ui.control.Button

        UserColumns string = string.empty   % user-chosen Units columns ([] = preset)
        ColPicker      matlab.ui.Figure       % the column-chooser window

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
        % Filter-bar columns (display labels come from filterLabel).
        FilterCols = ["region", "task", "bombcell_label", "c4_celltype"]
        % Preferred column order for the Units view (intersected with what's there).
        UnitsPreset = ["recording_tag", "unit", "cluster_id", "region", "task", ...
            "bombcell_label", "curation_label", "curation_note", "single_unit", ...
            "snr", "isi_viol_2p0_pct", "presence_ratio", "n_spikes", ...
            "purkinje_flag", "cs_rate_hz"]
        % Review states (in order) and their row background colours.
        Statuses = ["unverified", "WIP", "issue", "discard", "done"]
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
            if tool.appValid()
                tool.App.figureHandle().DeleteFcn = '';   % drop our close hook
            end
            if ~isempty(tool.ColPicker) && isvalid(tool.ColPicker)
                delete(tool.ColPicker);
            end
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

            barA = uigridlayout(controls, [1 9]);
            barA.Layout.Row = 1;
            barA.ColumnWidth = {120, "1x", 120, "1x", 90, 70, 55, 80, 130};
            barA.Padding = [0 0 0 0];
            uibutton(barA, Text="Load dataset...", ...
                ButtonPushedFcn=@(~, ~) tool.loadDatasetDialog());
            tool.DatasetLabel = uilabel(barA, Text="(no dataset)", ...
                FontColor=[0.3 0.3 0.3]);
            uibutton(barA, Text="Sortings root...", ...
                ButtonPushedFcn=@(~, ~) tool.chooseRootDialog());
            tool.RootLabel = uilabel(barA, Text="", FontColor=[0.3 0.3 0.3]);
            uibutton(barA, Text="Columns...", ...
                ButtonPushedFcn=@(~, ~) tool.chooseColumnsDialog());
            uibutton(barA, Text="Refresh", ...
                Tooltip="Re-read saved sortings: update curation labels/notes " + ...
                "and add units you created", ...
                ButtonPushedFcn=@(~, ~) tool.refreshFromDisk());
            tool.LockBox = uicheckbox(barA, Text="Lock", ...
                Tooltip="Lock review state (no changes are saved)", ...
                ValueChangedFcn=@(s, ~) tool.setLocked(s.Value));
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

            right = uigridlayout(mid, [6 1]);
            right.RowHeight = {"fit", "1x", "fit", "fit", "fit", 30};
            right.RowSpacing = 4;
            right.Padding = [0 0 0 0];
            uilabel(right, Text="Selected entry", FontWeight="bold");
            tool.DetailArea = uitextarea(right, Editable="off", Value="");
            uilabel(right, Text="Set review status", FontWeight="bold");
            sb = uigridlayout(right, [1 numel(tool.Statuses)]);
            sb.Padding = [0 0 0 0];
            sb.ColumnSpacing = 3;
            tool.StatusButtons = matlab.ui.control.Button.empty;
            for k = 1:numel(tool.Statuses)
                st = tool.Statuses(k);
                b = uibutton(sb, Text=st, ...
                    BackgroundColor=tool.statusColor(st), ...
                    FontColor=tool.statusFontColor(st), ...
                    ButtonPushedFcn=@(~, ~) tool.setStatusOfSelected(st));
                tool.StatusButtons(k) = b;
            end
            noteRow = uigridlayout(right, [1 2]);
            noteRow.Padding = [0 0 0 0];
            noteRow.ColumnWidth = {"fit", "1x"};
            uilabel(noteRow, Text="note", HorizontalAlignment="right");
            tool.NoteField = uieditfield(noteRow, "text", ...
                ValueChangedFcn=@(s, ~) tool.setNoteOfSelected(s.Value));
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
            tool.augmentTable();   % add cluster_id + curation_label/note columns
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
            tool.loadStoredColumns();
            tool.StateMap = containers.Map(KeyType="char", ValueType="any");
            tool.loadSidecar();
            tool.populateFilters();
            tool.applyFilters();
        end

        function loadStoredColumns(tool)
            %LOADSTOREDCOLUMNS Restore this dataset's chosen Units columns, if any.
            tool.UserColumns = string.empty;
            key = tool.colsPrefName();
            if key ~= "" && ispref(tool.PrefGroup, key)
                stored = string(getpref(tool.PrefGroup, key));
                tool.UserColumns = intersect(stored, tool.Dataset.vars, "stable");
            end
        end

        function key = colsPrefName(tool)
            key = "";
            if isfield(tool.Dataset, "source")
                [~, base] = fileparts(string(tool.Dataset.source));
                key = "cols_" + string(matlab.lang.makeValidName(base));
            end
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
                vars = tool.unitsDisplayVars();
            end
        end

        function vars = unitsDisplayVars(tool)
            %UNITSDISPLAYVARS Columns shown in the Units view: the user's chosen
            %   set if any, else the built-in preset (both intersected with what
            %   the dataset actually has).
            if ~isempty(tool.UserColumns)
                vars = intersect(tool.UserColumns, tool.Dataset.vars, "stable");
            else
                vars = intersect(tool.UnitsPreset, tool.Dataset.vars, "stable");
            end
            if isempty(vars)
                vars = tool.Dataset.vars;
            end
        end

        function chooseColumnsDialog(tool)
            %CHOOSECOLUMNSDIALOG Pick which columns the Units view shows; the
            %   choice is remembered per dataset.
            if isempty(tool.RawTable)
                tool.setInfo("Load a dataset first.");
                return
            end
            if ~isempty(tool.ColPicker) && isvalid(tool.ColPicker)
                figure(tool.ColPicker);
                return
            end
            allVars = tool.Dataset.vars;
            current = intersect(tool.unitsDisplayVars(), allVars, "stable");
            f = uifigure(Name="Choose Units columns", Position=[260 220 300 420]);
            tool.ColPicker = f;
            gl = uigridlayout(f, [3 2]);
            gl.RowHeight = {"fit", "1x", 30};
            gl.ColumnWidth = {"1x", "1x"};
            hint = uilabel(gl, Text="Select the columns to show (Units view):");
            hint.Layout.Row = 1; hint.Layout.Column = [1 2];
            lb = uilistbox(gl, Items=cellstr(allVars), Multiselect="on", ...
                Value=cellstr(current));
            lb.Layout.Row = 2; lb.Layout.Column = [1 2];
            ab = uibutton(gl, Text="Apply", BackgroundColor=[0.85 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.applyColumns(lb));
            ab.Layout.Row = 3; ab.Layout.Column = 1;
            rb = uibutton(gl, Text="Reset to default", ...
                ButtonPushedFcn=@(~, ~) tool.resetColumns());
            rb.Layout.Row = 3; rb.Layout.Column = 2;
        end

        function applyColumns(tool, lb)
            sel = string(lb.Value);
            tool.UserColumns = intersect(tool.Dataset.vars, sel, "stable");
            key = tool.colsPrefName();
            if key ~= "" && ~isempty(tool.UserColumns)
                setpref(tool.PrefGroup, key, cellstr(tool.UserColumns));
            end
            tool.closeColPicker();
            tool.renderTable();
        end

        function resetColumns(tool)
            tool.UserColumns = string.empty;
            key = tool.colsPrefName();
            if key ~= "" && ispref(tool.PrefGroup, key)
                rmpref(tool.PrefGroup, key);
            end
            tool.closeColPicker();
            tool.renderTable();
        end

        function closeColPicker(tool)
            if ~isempty(tool.ColPicker) && isvalid(tool.ColPicker)
                delete(tool.ColPicker);
            end
        end

        function renderTable(tool)
            if isempty(tool.ViewTable)
                tool.Table.Data = table();
                tool.setInfo("No rows match the current filters.");
                return
            end
            vars = tool.displayVars();
            n = height(tool.ViewTable);
            tags = string(tool.ViewTable.recording_tag);

            % Per-row review status, and a marker for the row loaded in the app.
            status = strings(n, 1);
            for i = 1:n
                [lvl, tg, un] = tool.entryOfRow(i);
                status(i) = tool.getStatus(lvl, tg, un);
            end
            isLoaded = tool.LoadedTag ~= "" & tags == tool.LoadedTag;
            if tool.RowKind == "units" && tool.LoadedUnit ~= "" ...
                    && ismember("unit", string(tool.ViewTable.Properties.VariableNames))
                isLoaded = isLoaded & string(tool.ViewTable.unit) == tool.LoadedUnit;
            end
            marker = strings(n, 1);
            marker(isLoaded) = "▶";

            data = [table(marker, status, VariableNames=["open", "review"]), ...
                tool.ViewTable(:, vars)];
            tool.Table.Data = data;
            tool.Table.ColumnName = [{''}, {'review'}, cellstr(vars(:)')];
            tool.Table.ColumnWidth = [{26}, {78}, repmat({'auto'}, 1, numel(vars))];

            removeStyle(tool.Table);
            % Colour each row by its review status.
            for k = 1:numel(tool.Statuses)
                st = tool.Statuses(k);
                rows = find(status == st);
                if isempty(rows)
                    continue
                end
                addStyle(tool.Table, uistyle(BackgroundColor=tool.statusColor(st)), ...
                    "row", rows);
                fc = tool.statusFontColor(st);
                if ~isequal(fc, [0 0 0])
                    addStyle(tool.Table, uistyle(FontColor=fc), "row", rows);
                end
            end
            % Grey out rows whose sorting folder (and online fallback) is absent.
            missing = arrayfun(@(t) ~tool.hasSorting(t), tags);
            if any(missing)
                addStyle(tool.Table, uistyle(FontColor=[0.6 0.6 0.6]), ...
                    "row", find(missing));
            end
            % Blue marker for the row currently loaded in the app.
            lr = find(isLoaded);
            if ~isempty(lr)
                addStyle(tool.Table, uistyle(FontColor=[0 0.35 0.9], ...
                    FontWeight="bold"), "cell", [lr(:), ones(numel(lr), 1)]);
            end
            lockMsg = "";
            if tool.Locked
                lockMsg = "  [LOCKED]";
            end
            tool.setInfo(sprintf("%d rows shown (%s).%s", n, tool.RowKind, lockMsg));
        end

        function onRowSelected(tool, ~)
            row = tool.selectedRow();
            if isempty(row)
                tool.DetailArea.Value = "";
                tool.NoteField.Value = "";
                return
            end
            [lvl, tg, un] = tool.entryOfRow(row);
            v = tool.getEntry(tool.statusKey(lvl, tg, un));
            head = "review: " + v.status;
            if v.modified ~= ""
                head = head + "   (modified " + v.modified + ")";
            end
            vars = string(tool.ViewTable.Properties.VariableNames);
            lines = strings(numel(vars) + 2, 1);
            lines(1) = head;
            lines(2) = "";
            for i = 1:numel(vars)
                lines(i + 2) = vars(i) + ": " + string(tool.ViewTable.(vars(i))(row));
            end
            tool.DetailArea.Value = cellstr(lines);
            tool.NoteField.Value = char(v.note);
            tool.OpenButton.Enable = ...
                matlab.lang.OnOffSwitchState(tool.hasSorting(tg));
        end
    end

    % ----------------------------------------------- curation / new units
    methods (Access = private)
        function augmentTable(tool)
            %AUGMENTTABLE Add cluster_id + curation_label/note columns and fill
            %   them from the saved sortings (existing units only; new units are
            %   added on demand by refreshFromDisk).
            t = tool.RawTable;
            n = height(t);
            if ~ismember("cluster_id", string(t.Properties.VariableNames))
                t.cluster_id = nan(n, 1);
            end
            t.curation_label = strings(n, 1);
            t.curation_note = strings(n, 1);
            tool.RawTable = t;
            tool.Dataset.vars = string(t.Properties.VariableNames);
            tool.buildClusterIds();   % cheap: one ALL_UNITS.tsv read
            % Curation labels/notes + new units are filled on demand (Refresh),
            % to keep opening the browser fast.
        end

        function m = finalIdMap(tool)
            %FINALIDMAP (recording_tag|merged_label) -> final cluster id, from
            %   the aggregate ALL_UNITS.tsv (one read).
            m = containers.Map(KeyType="char", ValueType="double");
            p = fullfile(tool.SortingsRoot, "ALL_UNITS.tsv");
            if ~isfile(p)
                return
            end
            try
                a = readtable(p, FileType="text", Delimiter="\t", ...
                    VariableNamingRule="preserve", TextType="string");
            catch
                return
            end
            need = ["recording_tag", "merged_label", "final_unit_id"];
            if ~all(ismember(need, string(a.Properties.VariableNames)))
                return
            end
            for i = 1:height(a)
                m(char(a.recording_tag(i) + "|" + a.merged_label(i))) = ...
                    double(a.final_unit_id(i));
            end
        end

        function buildClusterIds(tool)
            if ~tool.Dataset.hasRecordingTag ...
                    || ~ismember("unit", tool.Dataset.vars)
                return
            end
            m = tool.finalIdMap();
            if m.Count == 0
                return
            end
            t = tool.RawTable;
            for i = 1:height(t)
                key = char(string(t.recording_tag(i)) + "|" + string(t.unit(i)));
                if isKey(m, key)
                    t.cluster_id(i) = m(key);
                end
            end
            tool.RawTable = t;
        end

        function scanCuration(tool, addNew)
            %SCANCURATION Pull per-cluster Label/Note from each recording's
            %   cluster_curation.csv onto its unit rows; when addNew, also append
            %   rows for labelled clusters that are not dataset units.
            if ~tool.Dataset.hasRecordingTag
                return
            end
            t = tool.RawTable;
            tags = unique(string(t.recording_tag), "stable");
            added = {};
            nLabelled = 0;
            for tg = reshape(tags, 1, [])
                cc = fullfile(tool.SortingsRoot, tg, "phy_postmerge", ...
                    "cluster_curation.csv");
                if ~isfile(cc)
                    continue
                end
                try
                    C = readtable(cc, TextType="string");
                catch
                    continue
                end
                cvars = string(C.Properties.VariableNames);
                if ~ismember("ClusterID", cvars)
                    continue
                end
                rowsOfTag = find(string(t.recording_tag) == tg);
                idsOfTag = t.cluster_id(rowsOfTag);
                for j = 1:height(C)
                    cid = double(C.ClusterID(j));
                    lbl = colOrEmpty(C, "Label", j);
                    nt = colOrEmpty(C, "Note", j);
                    hit = rowsOfTag(idsOfTag == cid);
                    if ~isempty(hit)
                        t.curation_label(hit) = lbl;
                        t.curation_note(hit) = nt;
                        if lbl ~= ""
                            nLabelled = nLabelled + 1;
                        end
                    elseif addNew && lbl ~= ""
                        added{end + 1} = tool.makeNewUnitRow(t, tg, cid, lbl, nt); %#ok<AGROW>
                    end
                end
            end
            if ~isempty(added)
                t = [t; vertcat(added{:})];
            end
            tool.RawTable = t;
            tool.Dataset.vars = string(t.Properties.VariableNames);
            if addNew
                tool.setInfo(sprintf("Refreshed from sortings: %d labelled " + ...
                    "unit(s), %d added unit(s).", nLabelled, numel(added)));
            end
        end

        function r = makeNewUnitRow(~, t, tg, cid, lbl, nt)
            r = t(1, :);
            for v = string(r.Properties.VariableNames)
                col = r.(v);
                if isnumeric(col)
                    r.(v)(1) = NaN;
                elseif isstring(col)
                    r.(v)(1) = "";
                elseif islogical(col)
                    r.(v)(1) = false;
                else
                    try %#ok<TRYNC>
                        r.(v)(1) = missing;
                    end
                end
            end
            r.recording_tag(1) = tg;
            tt = t(string(t.recording_tag) == tg, :);
            if ~isempty(tt)
                if ismember("region", string(t.Properties.VariableNames))
                    r.region(1) = tt.region(1);
                end
                if ismember("task", string(t.Properties.VariableNames))
                    r.task(1) = tt.task(1);
                end
            end
            r.unit(1) = "c" + string(cid);
            r.cluster_id(1) = cid;
            r.curation_label(1) = lbl;
            r.curation_note(1) = nt;
        end

        function refreshFromDisk(tool)
            if isempty(tool.RawTable)
                tool.setInfo("Load a dataset first.");
                return
            end
            t = tool.RawTable;
            t.curation_label(:) = "";   % re-derive from disk (source of truth)
            t.curation_note(:) = "";
            tool.RawTable = t;
            tool.scanCuration(true);
            tool.populateFilters();
            tool.applyFilters();
        end
    end

    % ------------------------------------------------------- review status
    methods (Access = private)
        function [lvl, tg, un] = entryOfRow(tool, i)
            tg = string(tool.ViewTable.recording_tag(i));
            if tool.RowKind == "recordings"
                lvl = "recording";
                un = "";
            else
                lvl = "unit";
                if ismember("unit", string(tool.ViewTable.Properties.VariableNames))
                    un = string(tool.ViewTable.unit(i));
                else
                    un = "";
                end
            end
        end

        function k = statusKey(~, lvl, tg, un)
            k = char(string(lvl) + "|" + string(tg) + "|" + string(un));
        end

        function v = getEntry(tool, key)
            v = struct(status="unverified", note="", modified="");
            if ~isempty(tool.StateMap) && isKey(tool.StateMap, key)
                s = tool.StateMap(key);
                if isfield(s, "status") && s.status ~= ""
                    v.status = string(s.status);
                end
                if isfield(s, "note")
                    v.note = string(s.note);
                end
                if isfield(s, "modified")
                    v.modified = string(s.modified);
                end
            end
        end

        function s = getStatus(tool, lvl, tg, un)
            s = tool.getEntry(tool.statusKey(lvl, tg, un)).status;
        end

        function c = statusColor(~, st)
            switch string(st)
                case "WIP";     c = [1.00 0.93 0.55];
                case "issue";   c = [1.00 0.72 0.72];
                case "discard"; c = [0.32 0.32 0.32];
                case "done";    c = [0.78 0.92 0.78];
                otherwise;      c = [0.93 0.93 0.93];   % unverified
            end
        end

        function c = statusFontColor(~, st)
            if string(st) == "discard"
                c = [1 1 1];
            else
                c = [0 0 0];
            end
        end

        function applyState(tool, lvl, tg, un, statusVal, noteVal)
            key = tool.statusKey(lvl, tg, un);
            v = tool.getEntry(key);
            if ~isempty(statusVal)
                v.status = string(statusVal);
            end
            if ~isempty(noteVal)
                v.note = string(noteVal);
            end
            v.modified = string(datetime("now", Format="yyyy-MM-dd HH:mm"));
            tool.StateMap(key) = v;
            tool.saveSidecar();
        end

        function setStatusOfSelected(tool, st)
            if tool.Locked
                tool.setInfo("Locked - untick Lock to change review status.");
                return
            end
            row = tool.selectedRow();
            if isempty(row)
                tool.setInfo("Select a row first.");
                return
            end
            [lvl, tg, un] = tool.entryOfRow(row);
            tool.applyState(lvl, tg, un, st, []);
            tool.renderTable();
            tool.Table.Selection = row;
        end

        function setNoteOfSelected(tool, txt)
            if tool.Locked
                return
            end
            row = tool.selectedRow();
            if isempty(row)
                return
            end
            [lvl, tg, un] = tool.entryOfRow(row);
            tool.applyState(lvl, tg, un, [], string(txt));
        end

        function markOpened(tool, lvl, tg, un)
            % Opening an entry marks it WIP unless it already has a verdict.
            if tool.Locked || tg == ""
                return
            end
            if tool.getStatus(lvl, tg, un) == "unverified"
                tool.applyState(lvl, tg, un, "WIP", []);
            end
        end

        function markModified(tool, lvl, tg, un)
            if tool.Locked || tg == ""
                return
            end
            st = tool.getStatus(lvl, tg, un);
            if st == "unverified"
                st = "WIP";
            end
            tool.applyState(lvl, tg, un, st, []);   % refresh the modified stamp
        end

        function onAppClosed(tool)
            if isempty(tool.UIFigure) || ~isvalid(tool.UIFigure)
                return
            end
            if tool.LoadedTag ~= ""
                tool.markModified(tool.LoadedLevel, tool.LoadedTag, tool.LoadedUnit);
            end
            tool.LoadedTag = "";
            tool.LoadedUnit = "";
            tool.LoadedLevel = "";
            if ~isempty(tool.ViewTable)
                tool.renderTable();
            end
        end

        function setLocked(tool, tf)
            tool.Locked = logical(tf);
            en = matlab.lang.OnOffSwitchState(~tool.Locked);
            for b = tool.StatusButtons
                b.Enable = en;
            end
            if ~isempty(tool.NoteField) && isvalid(tool.NoteField)
                tool.NoteField.Enable = en;
            end
            if ~isempty(tool.ViewTable)
                tool.renderTable();
            end
        end

        function p = sidecarPath(tool)
            p = "";
            if isfield(tool.Dataset, "source") && tool.Dataset.source ~= ""
                [d, base] = fileparts(string(tool.Dataset.source));
                p = fullfile(d, base + ".review.tsv");
            end
        end

        function loadSidecar(tool)
            p = tool.sidecarPath();
            if p == "" || ~isfile(p)
                return
            end
            try
                t = readtable(p, FileType="text", Delimiter="\t", ...
                    VariableNamingRule="preserve", TextType="string");
            catch
                return
            end
            vars = string(t.Properties.VariableNames);
            if ~all(ismember(["level", "recording_tag", "unit", "status"], vars))
                return
            end
            hasNote = ismember("note", vars);
            hasMod = ismember("modified", vars);
            for i = 1:height(t)
                key = tool.statusKey(t.level(i), t.recording_tag(i), ...
                    emptyIfMissing(t.unit(i)));
                v.status = string(t.status(i));
                v.note = "";
                v.modified = "";
                if hasNote
                    v.note = emptyIfMissing(t.note(i));
                end
                if hasMod
                    v.modified = emptyIfMissing(t.modified(i));
                end
                tool.StateMap(key) = v;
            end
        end

        function saveSidecar(tool)
            if tool.Locked
                return
            end
            p = tool.sidecarPath();
            if p == ""
                return
            end
            keys = tool.StateMap.keys;
            n = numel(keys);
            if n == 0
                if isfile(p)
                    try
                        delete(p);
                    catch
                    end
                end
                return
            end
            [level, tg, un, status, note, modified] = deal(strings(n, 1));
            for i = 1:n
                parts = split(string(keys{i}), "|");
                level(i) = parts(1);
                tg(i) = parts(2);
                if numel(parts) >= 3
                    un(i) = parts(3);
                end
                v = tool.StateMap(keys{i});
                status(i) = fieldOr(v, "status");
                note(i) = fieldOr(v, "note");
                modified(i) = fieldOr(v, "modified");
            end
            T = table(level, tg, un, status, note, modified, VariableNames=...
                ["level", "recording_tag", "unit", "status", "note", "modified"]);
            try
                writetable(T, p, FileType="text", Delimiter="\t");
            catch err
                tool.setInfo("Could not write review sidecar: " + err.message);
            end
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
                cid = NaN;
                if ismember("cluster_id", string(tool.ViewTable.Properties.VariableNames))
                    cid = double(tool.ViewTable.cluster_id(row));
                end
                if isnan(cid)
                    cid = tool.resolveClusterId(tag, label);
                end
                if ~isnan(cid)
                    tool.App.setSelectedClusterIds(cid);
                    tool.App.openCuration(cid);   % skips silently if no waveforms
                else
                    tool.setInfo("Loaded " + tag + "; unit '" + label + ...
                        "' not found in final_labels.tsv (showing default cluster).");
                end
            end

            % Record what was loaded (for status + highlight) and hook the close.
            tool.LoadedTag = tag;
            if tool.RowKind == "units" ...
                    && ismember("unit", string(tool.ViewTable.Properties.VariableNames))
                tool.LoadedUnit = string(tool.ViewTable.unit(row));
                tool.LoadedLevel = "unit";
            else
                tool.LoadedUnit = "";
                tool.LoadedLevel = "recording";
            end
            tool.markOpened(tool.LoadedLevel, tool.LoadedTag, tool.LoadedUnit);
            if tool.appValid()
                tool.App.figureHandle().DeleteFcn = @(~, ~) tool.onAppClosed();
            end
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
        bombcell_label="Classification", c4_celltype="Cell type", ...
        purkinje_flag="Purkinje");
    if isfield(map, col)
        s = map.(col);
    else
        s = string(col);
    end
end

function s = emptyIfMissing(x)
    if ismissing(x)
        s = "";
    else
        s = string(x);
    end
end

function s = fieldOr(v, name)
    if isfield(v, name) && ~ismissing(v.(name))
        s = string(v.(name));
    else
        s = "";
    end
end

function s = colOrEmpty(t, name, i)
    if ismember(name, string(t.Properties.VariableNames))
        s = emptyIfMissing(t.(name)(i));
    else
        s = "";
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
