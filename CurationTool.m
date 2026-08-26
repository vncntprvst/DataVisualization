classdef CurationTool < handle
    %CURATIONTOOL Split or merge clusters using ISI, PC features and waveforms.
    %   CurationTool(app, clusterIds) opens a companion window for the given
    %   SpikeVisualizationApp and one or more clusters, showing three linked
    %   panels, each overlaying every selected cluster in its own colour:
    %     * ISI histogram (log axis)
    %     * PC-feature scatter (joint PCA of the selected clusters' waveforms)
    %     * the waveforms
    %
    %   Spikes can be selected three ways, across all shown clusters, and the
    %   selection appears in every panel (black):
    %     * "Circle ISI bars"  - freehand over the bars near 0; selects the
    %       spikes taking part in those short intervals.
    %     * "Lasso PC features" - freehand in the scatter; refines or selects.
    %     * "Line select"      - draw a line across the waveform panel; selects
    %       every shown spike whose waveform crosses it.
    %
    %   Then "Split off selected" moves the selection to a new cluster, or
    %   "Merge shown" combines the shown clusters (handy after overlaying two
    %   clusters to confirm they are one unit).
    %
    %   CurationTool(app) uses the app's currently selected cluster(s).
    %
    %   See also SpikeVisualizationApp.

    properties (Access = private)
        App                    % parent SpikeVisualizationApp (handle)
        ClusterIds double      % clusters shown (one or more)

        GlobalIdx double       % global spike indices shown (time order)
        ClusterOf double       % cluster id per shown spike
        IsiNextMs double       % ISI to the next spike in the same cluster, ms
        Scores double          % nSpikes-by-3 joint PCA scores (displayed)
        ScoresFull double      % nSpikes-by-(up to 10) PCA scores (clustered)
        TsneXY double = []     % cached 2-D t-SNE embedding ([] until computed)
        WaveUV double          % nSpikes-by-nSamples waveforms in microvolts
        WaveTimeMs double      % 1-by-nSamples waveform sample times, ms

        SelMask logical        % current spike selection (over shown spikes)
        SelSource string = ""  % "isi" | "pc" | "line" | "found"
        SelectedRange double = []   % [lo hi] ms of the circled ISI bars
        LineSegs cell = {}     % placed waveform lines; spikes must cross ALL
        FoundLabels double = []      % k-means labels from Find cluster ([] = none)

        ClusterMethod string = "gmm"      % "kmeans" | "gmm" | "dbscan"
        ClusterK = "auto"                 % "auto" or a positive integer
        ClusterPCs double = 6             % how many PCs feed Find clusters
        SettingsMenus struct = struct()

        UIFigure               matlab.ui.Figure
        ISIAxes                matlab.ui.control.UIAxes
        PCAxes                 matlab.ui.control.UIAxes
        WaveAxes               matlab.ui.control.UIAxes
        PCDropdown             matlab.ui.control.DropDown
        MoveDropdown           matlab.ui.control.DropDown
        FoundDropdown          matlab.ui.control.DropDown
        InfoLabel              matlab.ui.control.Label
        WindowTool                          % TimeWindowTool companion
        RecoverTool                         % RecoverSpikesTool companion
    end

    properties (Constant, Access = private)
        PrefGroup = "SpikeVisualizationApp"
        MinMs = 0.2
        MaxMs = 1000
        NumBins = 64
        RefractoryMs = 2
        MaxWaveformsToPlot = 150
        SelectColor = [0 0 0]
        PCChoices = [3 6 10]     % selectable PC counts for Find clusters
        MaxTsneSpikes = 8000     % t-SNE embeds at most this many spikes
    end

    methods
        function tool = CurationTool(app, clusterIds)
            arguments
                app (1,1) SpikeVisualizationApp
                clusterIds double = []
            end
            tool.App = app;
            if isempty(clusterIds)
                clusterIds = app.selectedClusterIds();
            end
            tool.ClusterIds = unique(clusterIds(:))';
            tool.computeData();
            tool.buildUI();
            % Restore the time-window split if it was open at last close.
            if getpref(tool.PrefGroup, "curationWindowTool", false)
                tool.launchWindowTool();
            end
            tool.refreshAll();
        end

        function delete(tool)
            % Remember whether the time-window split was open, so the next
            % CurationTool can restore it.
            try
                setpref(tool.PrefGroup, "curationWindowTool", ...
                    ~isempty(tool.WindowTool) && isvalid(tool.WindowTool));
            catch
                % prefs are a convenience - never block closing
            end
            if ~isempty(tool.WindowTool) && isvalid(tool.WindowTool)
                delete(tool.WindowTool);
            end
            if ~isempty(tool.RecoverTool) && isvalid(tool.RecoverTool)
                delete(tool.RecoverTool);
            end
            if ~isempty(tool.UIFigure) && isvalid(tool.UIFigure)
                delete(tool.UIFigure);
            end
        end
    end

    % ------------------------------------------- public (scriptable) selection
    methods
        function selectIsiRange(tool, rangeMs)
            %SELECTISIRANGE Select spikes taking part in ISIs within [lo hi] ms.
            arguments
                tool
                rangeMs (1, 2) double
            end
            tool.SelectedRange = rangeMs;
            tool.LineSegs = {};
            mask = false(numel(tool.GlobalIdx), 1);
            for cid = tool.ClusterIds
                ci = find(tool.ClusterOf == cid);           % time-ordered
                inR = tool.IsiNextMs(ci) >= rangeMs(1) & tool.IsiNextMs(ci) <= rangeMs(2);
                mask(ci(inR)) = true;                         % earlier spike
                mask(ci([false; inR(1:end-1)])) = true;       % later spike
            end
            tool.SelMask = mask;
            tool.SelSource = "isi";
            tool.setInfo(sprintf("ISI in [%.2f, %.2f] ms -> %d spikes selected.", ...
                rangeMs(1), rangeMs(2), sum(mask)));
            tool.refreshAll();
        end

        function selectPCLasso(tool, polyXY)
            %SELECTPCLASSO Select/refine spikes inside a polygon in the PC axes.
            arguments
                tool
                polyXY (:, 2) double
            end
            tool.LineSegs = {};
            [x, y] = tool.viewCoords();
            inside = inpolygon(x, y, polyXY(:, 1), polyXY(:, 2));
            if any(tool.SelMask)
                tool.SelMask = inside & tool.SelMask;
            else
                tool.SelMask = inside;
            end
            tool.SelSource = "pc";
            tool.SelectedRange = [];
            tool.setInfo(sprintf("%d spikes selected in the scatter.", ...
                sum(tool.SelMask)));
            tool.refreshAll();
        end

        function selectByLine(tool, p1, p2)
            %SELECTBYLINE Add a waveform line; select spikes crossing ALL lines.
            %   Call repeatedly to place several lines; a spike is selected only
            %   if its waveform crosses every one of them (a corridor).
            arguments
                tool
                p1 (1, 2) double
                p2 (1, 2) double
            end
            tool.LineSegs{end + 1} = [p1; p2];
            tool.SelSource = "line";
            tool.SelectedRange = [];
            tool.recomputeLineSelection();
            tool.refreshAll();
        end

        function clearLines(tool)
            %CLEARLINES Remove placed waveform lines and any found clustering.
            tool.LineSegs = {};
            tool.FoundLabels = [];
            tool.updateFoundDropdown();
            tool.SelMask = false(numel(tool.GlobalIdx), 1);
            tool.SelSource = "";
            tool.setInfo("Cleared lines and found clusters.");
            tool.refreshAll();
        end

        function newId = splitSelected(tool)
            %SPLITSELECTED Move the current selection into a new app cluster.
            newId = [];
            if ~any(tool.SelMask)
                tool.setInfo("Nothing selected.");
                return
            end
            newId = max(tool.App.Spikes.clusterIds) + 1;
            tool.App.applySplit(tool.GlobalIdx(tool.SelMask), newId);
            tool.ClusterIds = unique([tool.ClusterIds, newId]);
            tool.refreshAfterEdit();
            tool.setInfo(sprintf("Moved %d spikes to new cluster %d.", ...
                sum(tool.SelMask), newId));
        end

        function newId = moveSelected(tool)
            %MOVESELECTED Reassign the selected spikes to the destination cluster.
            %   Destination comes from the Move dropdown: an existing cluster, or
            %   "(new cluster)" to split them into a fresh one. This is how you
            %   merge just the selected spikes into another cluster.
            newId = [];
            if ~any(tool.SelMask)
                tool.setInfo("Nothing selected to move.");
                return
            end
            dest = tool.MoveDropdown.Value;   % [] => new cluster, else cluster id
            if isempty(dest)
                newId = max(tool.App.Spikes.clusterIds) + 1;
            else
                newId = dest;
            end
            n = sum(tool.SelMask);
            tool.App.applySplit(tool.GlobalIdx(tool.SelMask), newId);
            tool.ClusterIds = unique([tool.ClusterIds, newId]);
            tool.refreshAfterEdit();
            tool.setInfo(sprintf("Moved %d selected spikes to cluster %d.", n, newId));
        end

        function newId = mergeShown(tool)
            %MERGESHOWN Merge all the shown clusters into one.
            targets = unique(tool.ClusterIds(:));
            newId = tool.App.applyMerge(targets);
            if ~isempty(newId)
                tool.ClusterIds = newId;
                tool.refreshAfterEdit();
                tool.setInfo(sprintf("Merged %s into cluster %d.", ...
                    num2str(targets'), newId));
            end
        end

        function findClusters(tool)
            %FINDCLUSTERS Auto-split the shown spikes across PC features.
            %   Colours the scatter by the found clusters (overriding the source
            %   clusters); pick one in the Found dropdown to select/move it.
            X = tool.ScoresFull(:, 1:min(tool.ClusterPCs, size(tool.ScoresFull, 2)));
            if size(X, 1) < 10
                tool.setInfo("Too few spikes to cluster.");
                return
            end
            if tool.ClusterMethod == "dbscan"
                tool.setInfo("Finding clusters (dbscan)...");
                drawnow;
                tool.FoundLabels = tool.runDbscan(X);
            else
                k = tool.resolveK(X);
                tool.setInfo(sprintf("Finding %d clusters (%s)...", ...
                    k, tool.ClusterMethod));
                drawnow;
                tool.FoundLabels = tool.runClustering(X, k);
            end
            tool.SelMask = false(numel(tool.GlobalIdx), 1);
            tool.SelSource = "found";
            tool.SelectedRange = [];
            tool.LineSegs = {};
            tool.pickBestView(X);
            tool.updateFoundDropdown();
            tool.FoundDropdown.Value = tool.FoundDropdown.ItemsData{1};   % (all)
            nFound = numel(unique(tool.FoundLabels(tool.FoundLabels > 0)));
            msg = sprintf("Found %d clusters (%s, %d PCs).", nFound, ...
                tool.ClusterMethod, size(X, 2));
            nNoise = sum(tool.FoundLabels < 0);
            if nNoise > 0
                msg = msg + sprintf(" %d spikes flagged as noise.", nNoise);
            end
            tool.setInfo(msg + " Pick one in 'Found' to select it, then Move.");
            tool.refreshAll();
        end

        function refreshWindow(tool)
            %REFRESHWINDOW Recompute after the curation time window changed.
            tool.SelSource = "";
            tool.SelectedRange = [];
            tool.LineSegs = {};
            tool.FoundLabels = [];
            tool.computeData();
            tool.updateFoundDropdown();
            tool.refreshAll();
            tw = tool.App.timeWindowSec();
            tool.setInfo(sprintf("Time window %.1f-%.1f s: %d spikes shown.", ...
                tw(1), tw(2), numel(tool.GlobalIdx)));
        end

        function ids = shownClusters(tool)
            ids = tool.ClusterIds;
        end

        function g = globalSelection(tool)
            %GLOBALSELECTION App spike indices currently selected in this tool.
            g = tool.GlobalIdx(tool.SelMask);
        end

        function selectGlobal(tool, globalIdx)
            %SELECTGLOBAL Select the given app spike indices (from a linked tool).
            tool.SelMask = ismember(tool.GlobalIdx, globalIdx);
            tool.SelSource = "external";
            tool.setInfo(sprintf("%d spikes selected (from time-window lasso).", ...
                sum(tool.SelMask)));
            tool.refreshAll();
        end
    end

    methods (Access = private)
        function computeData(tool)
            tool.ClusterIds = reshape(tool.ClusterIds, 1, []);   % keep it a row
            s = tool.App.Spikes;
            inWin = ismember(s.clusters, tool.ClusterIds);
            tw = tool.App.timeWindowSec();
            if tw(1) > 0 || tw(2) < s.numSamples / s.samplingRate
                tSec = s.spikeTimes / s.samplingRate;
                inWin = inWin & tSec >= tw(1) & tSec <= tw(2);
            end
            idx = find(inWin);
            [~, order] = sort(s.spikeTimes(idx));
            tool.GlobalIdx = idx(order);
            tool.ClusterOf = s.clusters(tool.GlobalIdx);

            tool.IsiNextMs = nan(numel(tool.GlobalIdx), 1);
            for cid = tool.ClusterIds
                ci = find(tool.ClusterOf == cid);
                st = s.spikeTimes(tool.GlobalIdx(ci));
                if numel(st) > 1
                    tool.IsiNextMs(ci(1:end-1)) = diff(st) / s.samplingRate * 1000;
                end
            end

            uv = 1;
            if isfield(s, "uvPerADC")
                uv = s.uvPerADC;
            end
            tool.WaveUV = double(s.waveforms(tool.GlobalIdx, :)) * uv;
            w = s.waveformWindow;
            t = (-w(1):w(2)) / s.samplingRate * 1000;
            if numel(t) ~= size(tool.WaveUV, 2)
                t = 1:size(tool.WaveUV, 2);
            end
            tool.WaveTimeMs = t;

            if size(tool.WaveUV, 1) >= 3 && ~isempty(tool.WaveUV)
                nc = min([10, size(tool.WaveUV, 2), size(tool.WaveUV, 1) - 1]);
                [~, sc] = pca(tool.WaveUV, NumComponents=nc);
            else
                sc = zeros(numel(tool.GlobalIdx), 3);
            end
            if size(sc, 2) < 3
                sc(:, end+1:3) = 0;
            end
            tool.ScoresFull = sc;
            tool.Scores = sc(:, 1:3);
            tool.TsneXY = [];   % embedding is stale for the new spike set

            tool.SelMask = false(numel(tool.GlobalIdx), 1);
            tool.FoundLabels = [];
        end

        function k = resolveK(tool, X)
            if ~isequal(tool.ClusterK, "auto")
                k = double(tool.ClusterK);
                return
            end
            Xs = X(tool.subsample((1:size(X, 1))', 5000), :);   % keep it fast
            kMax = min(8, floor(size(Xs, 1) / 20));
            if kMax < 2
                k = 2;
                return
            end
            clustName = tool.ClusterMethod;
            if clustName == "gmm"
                clustName = "gmdistribution";
            end
            try
                % Calinski-Harabasz is O(n) (silhouette would be O(n^2)).
                ev = evalclusters(Xs, char(clustName), "CalinskiHarabasz", ...
                    KList=2:kMax);
                k = ev.OptimalK;
            catch
                k = 2;   % fall back if evaluation fails
            end
            if isnan(k) || k < 2
                k = 2;
            end
        end

        function labels = runClustering(tool, X, k)
            if tool.ClusterMethod == "gmm"
                try
                    gm = fitgmdist(X, k, RegularizationValue=1e-5, ...
                        Replicates=3, Options=statset(MaxIter=300));
                    labels = cluster(gm, X);
                    return
                catch
                    % fall through to k-means if the mixture fails to fit
                end
            end
            labels = kmeans(X, k, Replicates=3);
        end

        function labels = runDbscan(~, X)
            %RUNDBSCAN Density-based clustering (no k needed; -1 = noise).
            %   Epsilon comes from the standard kNN-distance heuristic: the
            %   90th percentile of each point's distance to its minpts-th
            %   neighbour, which sits near the elbow of the k-distance curve.
            minpts = max(10, 2 * size(X, 2));
            [~, d] = knnsearch(X, X, K=minpts);
            epsilon = quantile(d(:, end), 0.90);
            if ~(epsilon > 0)
                epsilon = max(eps, std(X(:)) / 10);   % degenerate features
            end
            labels = dbscan(X, epsilon, minpts);
        end

        function pickBestView(tool, X)
            %PICKBESTVIEW Show the 2-D PC projection that best separates the
            %   found clusters (auto-switches the PC dropdown if a better view
            %   than the current one exists).
            if string(tool.PCDropdown.Value) == "t-SNE"
                return   % keep the embedding view; found colours show there
            end
            pairs = [1 2; 1 3; 2 3];
            names = ["PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3"];
            sub = tool.subsample((1:size(X, 1))', 3000);   % silhouette is O(n^2)
            Xs = X(sub, :);
            labs = tool.FoundLabels(sub);
            score = -inf(1, size(pairs, 1));
            for p = 1:size(pairs, 1)
                if max(pairs(p, :)) > size(X, 2)
                    continue
                end
                try
                    score(p) = mean(silhouette(Xs(:, pairs(p, :)), labs));
                catch
                    score(p) = -inf;
                end
            end
            [~, best] = max(score);
            tool.PCDropdown.Value = names(best);
        end

        function selectFound(tool, g)
            if isempty(tool.FoundLabels)
                return
            end
            if isempty(g)
                tool.SelMask = false(numel(tool.GlobalIdx), 1);   % (all found)
            else
                tool.SelMask = tool.FoundLabels == g;
            end
            tool.SelSource = "found";
            tool.refreshAll();
        end

        function recomputeLineSelection(tool)
            if isempty(tool.LineSegs)
                tool.SelMask = false(numel(tool.GlobalIdx), 1);
            else
                mask = true(numel(tool.GlobalIdx), 1);
                for k = 1:numel(tool.LineSegs)
                    seg = tool.LineSegs{k};
                    mask = mask & spikesCrossingLine(tool.WaveUV, ...
                        tool.WaveTimeMs, seg(1, :), seg(2, :));
                end
                tool.SelMask = mask;
            end
            tool.setInfo(sprintf("%d line(s) placed; %d spikes cross all of them.", ...
                numel(tool.LineSegs), sum(tool.SelMask)));
        end

        function refreshAfterEdit(tool)
            tool.ClusterIds = intersect(tool.ClusterIds, tool.App.Spikes.clusterIds);
            tool.SelSource = "";
            tool.SelectedRange = [];
            tool.LineSegs = {};
            tool.FoundLabels = [];
            tool.computeData();
            tool.updateMoveDropdown();
            tool.updateFoundDropdown();
            tool.UIFigure.Name = tool.windowName();
            tool.refreshAll();
        end

        function refreshAll(tool)
            tool.refreshISI();
            tool.refreshPC();
            tool.refreshWave();
            if ~isempty(tool.WindowTool) && isvalid(tool.WindowTool)
                tool.WindowTool.setCurationSelection(tool.globalSelection());
            end
        end

        function name = windowName(tool)
            name = "Curation - clusters " + join(string(tool.ClusterIds), ", ");
        end

        function buildUI(tool)
            tool.UIFigure = uifigure(Name=tool.windowName(), ...
                Position=[120 110 1240 680], ...
                CloseRequestFcn=@(~, ~) delete(tool));   % close child tools too
            tool.buildSettingsMenu();
            outer = uigridlayout(tool.UIFigure, [4 1]);
            outer.RowHeight = {24, "1x", 78, 18};
            outer.Padding = [6 6 6 6];
            outer.RowSpacing = 4;

            % View controls above the panels (PC projection is over the scatter).
            head = uigridlayout(outer, [1 3]);
            head.Layout.Row = 1;
            head.Padding = [0 0 0 0];
            uilabel(head, Text="");
            pcCell = uigridlayout(head, [1 2]);
            pcCell.Layout.Column = 2;
            pcCell.ColumnWidth = {"fit", "1x"};
            pcCell.Padding = [0 0 0 0];
            uilabel(pcCell, Text="View:", HorizontalAlignment="right");
            tool.PCDropdown = uidropdown(pcCell, ...
                Items=["PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3", "t-SNE"], ...
                ValueChangedFcn=@(~, ~) tool.refreshPC());
            uilabel(head, Text="");

            axesRow = uigridlayout(outer, [1 3]);
            axesRow.Layout.Row = 2;
            axesRow.Padding = [0 0 0 0];
            tool.ISIAxes = uiaxes(axesRow);
            tool.PCAxes = uiaxes(axesRow);
            tool.WaveAxes = uiaxes(axesRow);

            % Button groups by category.
            groups = uigridlayout(outer, [1 3]);
            groups.Layout.Row = 3;
            groups.Padding = [0 0 0 0];
            groups.ColumnWidth = {"1.3x", "0.85x", "1.15x"};

            sel = uipanel(groups, Title="Selection");
            selg = uigridlayout(sel, [1 4]);
            selg.Padding = [4 2 4 2];
            uibutton(selg, Text="Circle ISI", ...
                ButtonPushedFcn=@(~, ~) tool.circleISI());
            uibutton(selg, Text="Lasso PC", ...
                ButtonPushedFcn=@(~, ~) tool.lassoPC());
            uibutton(selg, Text="Line select", ...
                ButtonPushedFcn=@(~, ~) tool.lineSelect());
            uibutton(selg, Text="Clear", ...
                ButtonPushedFcn=@(~, ~) tool.clearLines());

            clu = uipanel(groups, Title="Clustering");
            clug = uigridlayout(clu, [1 2]);
            clug.Padding = [4 2 4 2];
            uibutton(clug, Text="Find clusters", BackgroundColor=[0.9 0.9 1], ...
                Tooltip="Right-click to switch method (k-means / Gaussian mixture)", ...
                ContextMenu=tool.SettingsMenus.methodCtx, ...
                ButtonPushedFcn=@(~, ~) tool.findClusters());
            tool.FoundDropdown = uidropdown(clug, Items={'(select found)'}, ...
                ValueChangedFcn=@(s, ~) tool.selectFound(s.Value));

            act = uipanel(groups, Title="Actions");
            actg = uigridlayout(act, [1 3]);
            actg.Padding = [4 2 4 2];
            tool.MoveDropdown = uidropdown(actg, Items={'(new cluster)'});
            uibutton(actg, Text="Move selected", BackgroundColor=[1 0.85 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.moveSelected());
            uibutton(actg, Text="Merge shown", BackgroundColor=[0.85 1 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.mergeShown());

            tool.InfoLabel = uilabel(outer, Text="", FontColor=[0.2 0.2 0.2]);
            tool.InfoLabel.Layout.Row = 4;

            tool.updateMoveDropdown();
            tool.updateFoundDropdown();
            tool.setInfo("Select spikes (circle ISI / lasso PC / line), or " + ...
                "Find clusters to auto-split, then Move or Merge.");
        end

        function buildSettingsMenu(tool)
            m = uimenu(tool.UIFigure, Text="Settings");
            method = uimenu(m, Text="Find-cluster method");
            tool.SettingsMenus.kmeans = uimenu(method, Text="k-means", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("kmeans"));
            tool.SettingsMenus.gmm = uimenu(method, Text="Gaussian mixture", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("gmm"));
            tool.SettingsMenus.dbscan = uimenu(method, Text="DBSCAN (density)", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("dbscan"));
            % Same choices on right-click of the "Find clusters" button.
            cm = uicontextmenu(tool.UIFigure);
            tool.SettingsMenus.kmeansCtx = uimenu(cm, Text="k-means", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("kmeans"));
            tool.SettingsMenus.gmmCtx = uimenu(cm, Text="Gaussian mixture", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("gmm"));
            tool.SettingsMenus.dbscanCtx = uimenu(cm, Text="DBSCAN (density)", ...
                MenuSelectedFcn=@(~, ~) tool.setMethod("dbscan"));
            tool.SettingsMenus.methodCtx = cm;
            feat = uimenu(m, Text="Cluster on");
            tool.SettingsMenus.pcs = gobjects(1, numel(tool.PCChoices));
            for i = 1:numel(tool.PCChoices)
                c = tool.PCChoices(i);
                tool.SettingsMenus.pcs(i) = uimenu(feat, ...
                    Text=string(c) + " PCs", ...
                    MenuSelectedFcn=@(~, ~) tool.setClusterPCs(c));
            end
            tool.setMethod(tool.ClusterMethod);       % sync the check marks
            tool.setClusterPCs(tool.ClusterPCs);
            num = uimenu(m, Text="Number of clusters");
            tool.SettingsMenus.auto = uimenu(num, Text="Auto (estimate)", ...
                Checked="on", MenuSelectedFcn=@(~, ~) tool.setClusterK("auto"));
            tool.SettingsMenus.fixed = gobjects(1, 7);
            for kk = 2:8
                tool.SettingsMenus.fixed(kk-1) = uimenu(num, Text=string(kk), ...
                    MenuSelectedFcn=@(~, ~) tool.setClusterK(kk));
            end

            tools = uimenu(tool.UIFigure, Text="Tools");
            uimenu(tools, Text="Time-window split...", ...
                MenuSelectedFcn=@(~, ~) tool.launchWindowTool());
            uimenu(tools, Text="Recover missing spikes...", ...
                MenuSelectedFcn=@(~, ~) tool.launchRecoverTool());
        end

        function launchWindowTool(tool)
            if ~isempty(tool.WindowTool) && isvalid(tool.WindowTool)
                delete(tool.WindowTool);
            end
            tool.WindowTool = TimeWindowTool(tool.App, tool);
            tool.arrangeWithWindowTool();
        end

        function arrangeWithWindowTool(tool)
            %ARRANGEWITHWINDOWTOOL Tile this window and the time-window split
            %   side by side (Curation left, split right) on the monitor the
            %   Curation window is on, so neither covers the other.
            cf = tool.UIFigure;
            wf = tool.WindowTool.figureHandle();
            scr = tool.monitorOf(cf);
            gap = 10;                        % between/around the windows
            topPad = 90;                     % OS title bar above the figure
            botPad = 60;                     % taskbar
            availW = scr(3) - 3 * gap;
            availH = scr(4) - topPad - botPad;
            curW = min(cf.Position(3), round(availW * 0.55));
            wtW = min(wf.Position(3), availW - curW);
            yTop = scr(2) + botPad + availH;
            curH = min(cf.Position(4), availH);
            wtH = min(wf.Position(4), availH);
            cf.Position = [scr(1) + gap, yTop - curH, curW, curH];
            wf.Position = [scr(1) + 2 * gap + curW, yTop - wtH, wtW, wtH];
        end

        function scr = monitorOf(~, fig)
            %MONITOROF The monitor rect [x y w h] containing the figure centre.
            mons = groot().MonitorPositions;
            scr = mons(1, :);
            c = fig.Position(1:2) + fig.Position(3:4) / 2;
            for i = 1:size(mons, 1)
                p = mons(i, :);
                if c(1) >= p(1) && c(1) < p(1) + p(3) ...
                        && c(2) >= p(2) && c(2) < p(2) + p(4)
                    scr = p;
                    break
                end
            end
        end

        function launchRecoverTool(tool)
            if ~isempty(tool.RecoverTool) && isvalid(tool.RecoverTool)
                delete(tool.RecoverTool);
            end
            tool.RecoverTool = RecoverSpikesTool(tool.App, tool);
        end

        function setMethod(tool, method)
            tool.ClusterMethod = method;
            for nm = ["kmeans", "gmm", "dbscan"]
                on = matlab.lang.OnOffSwitchState(nm == method);
                tool.SettingsMenus.(nm).Checked = on;
                tool.SettingsMenus.(nm + "Ctx").Checked = on;
            end
        end

        function setClusterPCs(tool, n)
            tool.ClusterPCs = n;
            for i = 1:numel(tool.PCChoices)
                tool.SettingsMenus.pcs(i).Checked = ...
                    matlab.lang.OnOffSwitchState(tool.PCChoices(i) == n);
            end
        end

        function setClusterK(tool, k)
            tool.ClusterK = k;
            tool.SettingsMenus.auto.Checked = ...
                matlab.lang.OnOffSwitchState(isequal(k, "auto"));
            for kk = 2:8
                tool.SettingsMenus.fixed(kk-1).Checked = ...
                    matlab.lang.OnOffSwitchState(isequal(k, kk));
            end
        end

        function updateMoveDropdown(tool)
            ids = tool.App.Spikes.clusterIds;
            tool.MoveDropdown.Items = ["(new cluster)"; "-> cluster " + string(ids(:))];
            tool.MoveDropdown.ItemsData = [{[]}; num2cell(ids(:))];
        end

        function updateFoundDropdown(tool)
            if isempty(tool.FoundLabels)
                tool.FoundDropdown.Items = {'(select found)'};
                tool.FoundDropdown.ItemsData = [];
                return
            end
            g = unique(tool.FoundLabels);
            names = "found " + string(g(:));
            names(g < 0) = "noise";   % dbscan outliers
            tool.FoundDropdown.Items = ["(all found)"; names];
            tool.FoundDropdown.ItemsData = [{[]}; num2cell(g(:))];
        end

        function circleISI(tool)
            roi = drawfreehand(tool.ISIAxes, Color=[0.2 0.2 0.2]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 3
                return
            end
            xRange = [min(roi.Position(:, 1)), max(roi.Position(:, 1))];
            delete(roi);
            tool.selectIsiRange(xRange);
        end

        function lassoPC(tool)
            roi = drawfreehand(tool.PCAxes, Color=[0.2 0.2 0.2]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 3
                return
            end
            poly = roi.Position;
            delete(roi);
            tool.selectPCLasso(poly);
        end

        function lineSelect(tool)
            roi = drawline(tool.WaveAxes, Color=[0.2 0.2 0.2]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 2
                return
            end
            p = roi.Position;
            delete(roi);
            tool.selectByLine(p(1, :), p(2, :));
        end

        function refreshISI(tool)
            ax = tool.ISIAxes;
            cla(ax);
            hold(ax, "on");
            edges = logspace(log10(tool.MinMs), log10(tool.MaxMs), tool.NumBins);
            for cid = tool.ClusterIds
                isi = tool.IsiNextMs(tool.ClusterOf == cid);
                isi = isi(~isnan(isi));
                histogram(ax, isi, edges, FaceColor=tool.App.clusterColor(cid), ...
                    FaceAlpha=0.45, EdgeColor="none");
            end
            selIsi = tool.selectionIsi();
            if ~isempty(selIsi)
                histogram(ax, selIsi, edges, FaceColor=tool.SelectColor, ...
                    FaceAlpha=0.6, EdgeColor="none");
            end
            xline(ax, tool.RefractoryMs, "--", Color=[0.6 0 0]);
            set(ax, XScale="log");
            hold(ax, "off");
            allIsi = tool.IsiNextMs(~isnan(tool.IsiNextMs));
            pct = 100 * sum(allIsi < tool.RefractoryMs) / max(1, numel(allIsi));
            title(ax, sprintf("ISI  (%.2f%% < %g ms)", pct, tool.RefractoryMs));
            xlabel(ax, "ISI (ms, log)");
            ylabel(ax, "Count");
            xlim(ax, [tool.MinMs tool.MaxMs]);
        end

        function selIsi = selectionIsi(tool)
            if tool.SelSource == "isi" && ~isempty(tool.SelectedRange)
                isi = tool.IsiNextMs;
                selIsi = isi(isi >= tool.SelectedRange(1) & ...
                    isi <= tool.SelectedRange(2));
            elseif any(tool.SelMask)
                selIsi = tool.IsiNextMs(tool.SelMask & ~isnan(tool.IsiNextMs));
            else
                selIsi = [];
            end
        end

        function [masks, colors, names] = colorGroups(tool)
            %COLORGROUPS Groups for colouring: found clusters if a Find has run,
            %   otherwise the source clusters.
            if ~isempty(tool.FoundLabels)
                g = unique(tool.FoundLabels);
                masks = arrayfun(@(v) tool.FoundLabels == v, g, ...
                    UniformOutput=false);
                colors = foundPalette(numel(g));
                colors(g < 0, :) = 0.55;   % grey for dbscan noise
                names = "found " + string(g(:));
                names(g < 0) = "noise";
            else
                cids = tool.ClusterIds;
                masks = arrayfun(@(c) tool.ClusterOf == c, cids, ...
                    UniformOutput=false);
                colors = zeros(numel(cids), 3);
                for k = 1:numel(cids)
                    colors(k, :) = tool.App.clusterColor(cids(k));
                end
                names = "cluster " + string(cids(:));
            end
        end

        function refreshPC(tool)
            ax = tool.PCAxes;
            [x, y, xName, yName] = tool.viewCoords();
            cla(ax);
            hold(ax, "on");
            [masks, colors, names] = tool.colorGroups();
            handles = gobjects(1, numel(masks));
            for k = 1:numel(masks)
                handles(k) = scatter(ax, x(masks{k}), y(masks{k}), 6, ...
                    colors(k, :), "filled", MarkerFaceAlpha=0.35);
            end
            if any(tool.SelMask)
                scatter(ax, x(tool.SelMask), y(tool.SelMask), 14, ...
                    tool.SelectColor, "o", LineWidth=0.75);
            end
            hold(ax, "off");
            if numel(names) > 1
                legend(ax, handles, names, Location="best", Box="off");
            end
            if xName == "t-SNE 1"
                title(ax, sprintf("t-SNE embedding  (of %d PCs)", ...
                    size(tool.ScoresFull, 2)));
            else
                title(ax, "PC features  (" + string(tool.PCDropdown.Value) + ")");
            end
            xlabel(ax, xName);
            ylabel(ax, yName);
            axis(ax, "tight");
        end

        function [x, y, xName, yName] = viewCoords(tool)
            %VIEWCOORDS Scatter coordinates for the current view (PC pair or
            %   t-SNE). t-SNE coordinates are NaN for spikes left out of the
            %   embedding, which scatter and inpolygon both ignore.
            if string(tool.PCDropdown.Value) == "t-SNE"
                xy = tool.tsneCoords();
                x = xy(:, 1);
                y = xy(:, 2);
                xName = "t-SNE 1";
                yName = "t-SNE 2";
            else
                [ax1, ax2] = tool.pcColumns();
                x = tool.Scores(:, ax1);
                y = tool.Scores(:, ax2);
                parts = split(string(tool.PCDropdown.Value), " vs ");
                xName = parts(1);
                yName = parts(2);
            end
        end

        function xy = tsneCoords(tool)
            %TSNECOORDS Lazily computed 2-D t-SNE of the clustering PCs,
            %   cached until the spike set changes (computeData clears it).
            if ~isempty(tool.TsneXY)
                xy = tool.TsneXY;
                return
            end
            n = size(tool.ScoresFull, 1);
            xy = nan(n, 2);
            sub = tool.subsample((1:n)', tool.MaxTsneSpikes);
            if numel(sub) >= 12
                note = "";
                if numel(sub) < n
                    note = sprintf(" (%d of %d spikes)", numel(sub), n);
                end
                tool.setInfo(sprintf("Computing t-SNE embedding%s - this " + ...
                    "can take a while...", note));
                drawnow;
                try
                    xy(sub, :) = tsne(tool.ScoresFull(sub, :), ...
                        Perplexity=min(30, floor((numel(sub) - 1) / 3)));
                    tool.setInfo("t-SNE ready" + note + ". Lasso and Find " + ...
                        "clusters work in this view too.");
                catch err
                    tool.setInfo("t-SNE failed: " + err.message);
                end
            else
                tool.setInfo("Too few spikes for a t-SNE embedding.");
            end
            tool.TsneXY = xy;
        end

        function refreshWave(tool)
            ax = tool.WaveAxes;
            cla(ax);
            hold(ax, "on");
            t = tool.WaveTimeMs;
            [masks, colors] = tool.colorGroups();
            for k = 1:numel(masks)
                ci = tool.subsample(find(masks{k}), tool.MaxWaveformsToPlot);
                if ~isempty(ci)
                    plot(ax, t, tool.WaveUV(ci, :)', Color=[colors(k, :) 0.15]);
                end
            end
            sel = tool.subsample(find(tool.SelMask), tool.MaxWaveformsToPlot);
            if ~isempty(sel)
                plot(ax, t, tool.WaveUV(sel, :)', ...
                    Color=[tool.SelectColor 0.35], LineWidth=0.5);
            end
            for k = 1:numel(tool.LineSegs)   % show the placed selection lines
                seg = tool.LineSegs{k};
                plot(ax, seg(:, 1), seg(:, 2), "--", Color=[0.8 0.1 0.1], ...
                    LineWidth=1);
            end
            hold(ax, "off");
            title(ax, sprintf("Waveforms  (%d selected, %d line(s))", ...
                sum(tool.SelMask), numel(tool.LineSegs)));
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Voltage (\muV)");
            axis(ax, "tight");
        end

        function [ax1, ax2] = pcColumns(tool)
            switch tool.PCDropdown.Value
                case "PC1 vs PC3"
                    ax1 = 1; ax2 = 3;
                case "PC2 vs PC3"
                    ax1 = 2; ax2 = 3;
                otherwise
                    ax1 = 1; ax2 = 2;
            end
        end

        function setInfo(tool, msg)
            if ~isempty(tool.InfoLabel) && isvalid(tool.InfoLabel)
                tool.InfoLabel.Text = string(msg);
            end
        end
    end

    methods (Static, Access = private)
        function idx = subsample(idx, maxN)
            if numel(idx) > maxN
                idx = idx(round(linspace(1, numel(idx), maxN)));
            end
        end
    end
end

% =========================================================================
function pal = foundPalette(k)
    %FOUNDPALETTE k well-separated colours for found clusters (golden-angle).
    pal = zeros(k, 3);
    for i = 1:k
        pal(i, :) = hsv2rgb([mod((i - 1) * 0.6180339887, 1), 0.7, 0.85]);
    end
end

function mask = spikesCrossingLine(waveformsUV, timeMs, p1, p2)
    %SPIKESCROSSINGLINE Which waveforms intersect the segment p1->p2.
    %   waveformsUV is nSpikes-by-nSamples (microvolts) sampled at timeMs.
    %   p1, p2 are [timeMs voltage] endpoints. Returns an nSpikes logical.
    numSpikes = size(waveformsUV, 1);
    mask = false(numSpikes, 1);
    dx = p2(1) - p1(1);
    dy = p2(2) - p1(2);
    for k = 1:size(waveformsUV, 2) - 1
        ax = timeMs(k);
        bx = timeMs(k + 1);
        ay = waveformsUV(:, k);
        by = waveformsUV(:, k + 1);
        d1 = dx * (ay - p1(2)) - dy * (ax - p1(1));
        d2 = dx * (by - p1(2)) - dy * (bx - p1(1));
        ex = bx - ax;
        ey = by - ay;
        d3 = ex .* (p1(2) - ay) - ey .* (p1(1) - ax);
        d4 = ex .* (p2(2) - ay) - ey .* (p2(1) - ax);
        mask = mask | (d1 .* d2 < 0 & d3 .* d4 < 0);
    end
end
