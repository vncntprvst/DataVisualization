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
        Scores double          % nSpikes-by-3 joint PCA scores of the waveforms
        WaveUV double          % nSpikes-by-nSamples waveforms in microvolts
        WaveTimeMs double      % 1-by-nSamples waveform sample times, ms

        SelMask logical        % current spike selection (over shown spikes)
        SelSource string = ""  % "isi" | "pc" | "line"
        SelectedRange double = []   % [lo hi] ms of the circled ISI bars
        LineSegs cell = {}     % placed waveform lines; spikes must cross ALL

        UIFigure               matlab.ui.Figure
        ISIAxes                matlab.ui.control.UIAxes
        PCAxes                 matlab.ui.control.UIAxes
        WaveAxes               matlab.ui.control.UIAxes
        PCDropdown             matlab.ui.control.DropDown
        MoveDropdown           matlab.ui.control.DropDown
        MergeDropdown          matlab.ui.control.DropDown
        InfoLabel              matlab.ui.control.Label
    end

    properties (Constant, Access = private)
        MinMs = 0.2
        MaxMs = 1000
        NumBins = 64
        RefractoryMs = 2
        MaxWaveformsToPlot = 150
        SelectColor = [0 0 0]
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
            tool.refreshAll();
        end

        function delete(tool)
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
            [ax1, ax2] = tool.pcColumns();
            inside = inpolygon(tool.Scores(:, ax1), tool.Scores(:, ax2), ...
                polyXY(:, 1), polyXY(:, 2));
            if any(tool.SelMask)
                tool.SelMask = inside & tool.SelMask;
            else
                tool.SelMask = inside;
            end
            tool.SelSource = "pc";
            tool.SelectedRange = [];
            tool.setInfo(sprintf("%d spikes selected in PC space.", sum(tool.SelMask)));
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
            %CLEARLINES Remove all placed waveform lines and their selection.
            tool.LineSegs = {};
            if tool.SelSource == "line"
                tool.SelMask = false(numel(tool.GlobalIdx), 1);
                tool.SelSource = "";
            end
            tool.setInfo("Cleared waveform lines.");
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
            %MERGESHOWN Merge the shown clusters (plus the dropdown target).
            targets = tool.ClusterIds(:);
            if ~isempty(tool.MergeDropdown.ItemsData) ...
                    && ~isempty(tool.MergeDropdown.Value)
                targets = [targets; tool.MergeDropdown.Value];
            end
            targets = unique(targets);
            newId = tool.App.applyMerge(targets);
            if ~isempty(newId)
                tool.ClusterIds = newId;
                tool.refreshAfterEdit();
                tool.setInfo(sprintf("Merged %s into cluster %d.", ...
                    num2str(targets'), newId));
            end
        end
    end

    methods (Access = private)
        function computeData(tool)
            tool.ClusterIds = reshape(tool.ClusterIds, 1, []);   % keep it a row
            s = tool.App.Spikes;
            idx = find(ismember(s.clusters, tool.ClusterIds));
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
                [~, sc] = pca(tool.WaveUV, NumComponents=min(3, size(tool.WaveUV, 2)));
                if size(sc, 2) < 3
                    sc(:, end+1:3) = 0;
                end
            else
                sc = zeros(numel(tool.GlobalIdx), 3);
            end
            tool.Scores = sc;

            tool.SelMask = false(numel(tool.GlobalIdx), 1);
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
            tool.computeData();
            tool.updateMoveDropdown();
            tool.updateMergeDropdown();
            tool.UIFigure.Name = tool.windowName();
            tool.refreshAll();
        end

        function refreshAll(tool)
            tool.refreshISI();
            tool.refreshPC();
            tool.refreshWave();
        end

        function name = windowName(tool)
            name = "Curation - clusters " + join(string(tool.ClusterIds), ", ");
        end

        function buildUI(tool)
            tool.UIFigure = uifigure(Name=tool.windowName(), ...
                Position=[140 140 1200 620]);
            outer = uigridlayout(tool.UIFigure, [3 1]);
            outer.RowHeight = {"1x", 36, 20};
            outer.Padding = [6 6 6 6];

            axesRow = uigridlayout(outer, [1 3]);
            axesRow.Layout.Row = 1;
            axesRow.Padding = [0 0 0 0];
            tool.ISIAxes = uiaxes(axesRow);
            tool.PCAxes = uiaxes(axesRow);
            tool.WaveAxes = uiaxes(axesRow);

            btnRow = uigridlayout(outer, [1 9]);
            btnRow.Layout.Row = 2;
            btnRow.Padding = [0 0 0 0];
            uibutton(btnRow, Text="Circle ISI bars", ...
                ButtonPushedFcn=@(~, ~) tool.circleISI());
            uibutton(btnRow, Text="Lasso PC features", ...
                ButtonPushedFcn=@(~, ~) tool.lassoPC());
            uibutton(btnRow, Text="Line select (add)", ...
                ButtonPushedFcn=@(~, ~) tool.lineSelect());
            uibutton(btnRow, Text="Clear lines", ...
                ButtonPushedFcn=@(~, ~) tool.clearLines());
            tool.PCDropdown = uidropdown(btnRow, ...
                Items=["PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3"], ...
                ValueChangedFcn=@(~, ~) tool.refreshPC());
            tool.MoveDropdown = uidropdown(btnRow, Items={'(new cluster)'});
            uibutton(btnRow, Text="Move selected", ...
                BackgroundColor=[1 0.85 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.moveSelected());
            tool.MergeDropdown = uidropdown(btnRow, Items={'(none)'});
            uibutton(btnRow, Text="Merge shown", BackgroundColor=[0.85 1 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.mergeShown());

            tool.InfoLabel = uilabel(outer, Text="", FontColor=[0.2 0.2 0.2]);
            tool.InfoLabel.Layout.Row = 3;

            tool.updateMoveDropdown();
            tool.updateMergeDropdown();
            tool.setInfo("Select spikes by circling ISI bars, lassoing PC " + ...
                "features, or drawing a line across the waveforms.");
        end

        function updateMoveDropdown(tool)
            ids = tool.App.Spikes.clusterIds;
            tool.MoveDropdown.Items = ["(new cluster)"; "-> cluster " + string(ids(:))];
            tool.MoveDropdown.ItemsData = [{[]}; num2cell(ids(:))];
        end

        function updateMergeDropdown(tool)
            others = setdiff(tool.App.Spikes.clusterIds, tool.ClusterIds);
            if isempty(others)
                tool.MergeDropdown.Items = {'(also merge with...)'};
                tool.MergeDropdown.ItemsData = [];
            else
                tool.MergeDropdown.Items = ["(also merge with...)"; ...
                    "cluster " + string(others(:))];
                tool.MergeDropdown.ItemsData = [{[]}; num2cell(others(:))];
            end
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

        function refreshPC(tool)
            ax = tool.PCAxes;
            cla(ax);
            hold(ax, "on");
            [ax1, ax2] = tool.pcColumns();
            x = tool.Scores(:, ax1);
            y = tool.Scores(:, ax2);
            handles = gobjects(1, numel(tool.ClusterIds));
            for k = 1:numel(tool.ClusterIds)
                m = tool.ClusterOf == tool.ClusterIds(k);
                handles(k) = scatter(ax, x(m), y(m), 6, ...
                    tool.App.clusterColor(tool.ClusterIds(k)), "filled", ...
                    MarkerFaceAlpha=0.35);
            end
            if any(tool.SelMask)
                scatter(ax, x(tool.SelMask), y(tool.SelMask), 14, ...
                    tool.SelectColor, "o", LineWidth=0.75);
            end
            hold(ax, "off");
            if numel(tool.ClusterIds) > 1
                legend(ax, handles, "cluster " + string(tool.ClusterIds), ...
                    Location="best", Box="off");
            end
            labels = tool.PCDropdown.Value;
            title(ax, "PC features  (" + labels + ")");
            parts = split(labels, " vs ");
            xlabel(ax, parts(1));
            ylabel(ax, parts(2));
            axis(ax, "tight");
        end

        function refreshWave(tool)
            ax = tool.WaveAxes;
            cla(ax);
            hold(ax, "on");
            t = tool.WaveTimeMs;
            for cid = tool.ClusterIds
                ci = tool.subsample(find(tool.ClusterOf == cid), tool.MaxWaveformsToPlot);
                if ~isempty(ci)
                    plot(ax, t, tool.WaveUV(ci, :)', ...
                        Color=[tool.App.clusterColor(cid) 0.15]);
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
