classdef ISISplitTool < handle
    %ISISPLITTOOL Split refractory-violating spikes out of a cluster.
    %   ISISplitTool(app, clusterId) opens a companion window for the given
    %   SpikeVisualizationApp and cluster. It shows the cluster's inter-spike-
    %   interval (ISI) histogram next to a PC-feature scatter (temporal PCA of
    %   the waveforms, the same features Phy's FeatureView shows).
    %
    %   Workflow:
    %     1. "Circle ISI bars" - draw a freehand region over the bars near 0.
    %        Every spike that takes part in an ISI inside that range is
    %        highlighted (orange) in the PC-feature scatter.
    %     2. "Lasso PC features" - optionally draw a freehand region in the
    %        scatter to refine the exact spikes to remove (red).
    %     3. "Split off selected" - the highlighted spikes (the PC lasso if one
    %        was drawn, otherwise the ISI selection) are moved to a new cluster
    %        in the parent app, so you can judge and keep or discard them.
    %
    %   ISISplitTool(app) uses the app's currently selected cluster.
    %
    %   See also SpikeVisualizationApp.

    properties (Access = private)
        App                    % parent SpikeVisualizationApp (handle)
        ClusterId double

        GlobalIdx double       % global spike indices of this cluster (time order)
        IsiNextMs double       % ISI to the next spike, ms (NaN for the last)
        Scores double          % nSpikes-by-3 PCA scores of the waveforms

        IsiMask logical        % spikes taking part in the selected ISI bins
        PCMask logical         % spikes inside the PC-feature lasso

        UIFigure               matlab.ui.Figure
        ISIAxes                matlab.ui.control.UIAxes
        PCAxes                 matlab.ui.control.UIAxes
        PCDropdown             matlab.ui.control.DropDown
        InfoLabel              matlab.ui.control.Label
        Histogram              matlab.graphics.chart.primitive.Histogram
    end

    properties (Constant, Access = private)
        MaxMs = 50
        BinMs = 0.5
        RefractoryMs = 2
    end

    methods
        function tool = ISISplitTool(app, clusterId)
            arguments
                app (1,1) SpikeVisualizationApp
                clusterId double = []
            end
            tool.App = app;
            if isempty(clusterId)
                clusterId = app.currentCluster();
            end
            tool.ClusterId = clusterId;
            tool.computeData();
            tool.buildUI();
            tool.refreshISI();
            tool.refreshPC();
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
            %SELECTISIRANGE Highlight spikes whose ISI falls in [lo hi] ms.
            %   Flags both the earlier and the later spike of every ISI in range.
            arguments
                tool
                rangeMs (1, 2) double
            end
            inRange = tool.IsiNextMs >= rangeMs(1) & tool.IsiNextMs <= rangeMs(2);
            participants = false(numel(tool.GlobalIdx), 1);
            participants(inRange) = true;                     % earlier spike
            participants([false; inRange(1:end-1)]) = true;   % later spike
            tool.IsiMask = participants;
            tool.PCMask = false(numel(tool.GlobalIdx), 1);    % reset refinement
            tool.setInfo(sprintf("ISI in [%.2f, %.2f] ms -> %d spikes highlighted. " + ...
                "Lasso the PC scatter to refine, then split.", ...
                rangeMs(1), rangeMs(2), sum(tool.IsiMask)));
            tool.refreshISI();
            tool.refreshPC();
        end

        function selectPCLasso(tool, polyXY)
            %SELECTPCLASSO Select spikes inside a polygon in the current PC axes.
            %   polyXY is an N-by-2 list of [PCa PCb] vertices matching the view.
            arguments
                tool
                polyXY (:, 2) double
            end
            [ax1, ax2] = tool.pcColumns();
            inside = inpolygon(tool.Scores(:, ax1), tool.Scores(:, ax2), ...
                polyXY(:, 1), polyXY(:, 2));
            if any(tool.IsiMask)
                tool.PCMask = inside & tool.IsiMask;   % refine within ISI set
            else
                tool.PCMask = inside;
            end
            tool.setInfo(sprintf("%d spikes lassoed in PC space. " + ...
                "Click 'Split off selected' to move them.", sum(tool.PCMask)));
            tool.refreshPC();
        end

        function newId = splitSelected(tool)
            %SPLITSELECTED Move the active selection into a new app cluster.
            %   Uses the PC lasso if one was drawn, otherwise the ISI selection.
            newId = [];
            if any(tool.PCMask)
                mask = tool.PCMask;
            elseif any(tool.IsiMask)
                mask = tool.IsiMask;
            else
                tool.setInfo("Nothing selected. Circle ISI bars or lasso PC features first.");
                return
            end
            newId = max(tool.App.Spikes.clusterIds) + 1;
            globalToMove = tool.GlobalIdx(mask);
            tool.App.applySplit(globalToMove, newId);
            tool.computeData();       % source cluster is now smaller
            tool.refreshISI();
            tool.refreshPC();
            tool.setInfo(sprintf("Moved %d spikes to cluster %d. " + ...
                "Source cluster %d now has %d spikes.", numel(globalToMove), ...
                newId, tool.ClusterId, numel(tool.GlobalIdx)));
        end
    end

    methods (Access = private)
        function computeData(tool)
            %COMPUTEDATA Gather this cluster's spikes, ISIs and waveform PCA.
            s = tool.App.Spikes;
            idx = find(s.clusters == tool.ClusterId);
            [~, order] = sort(s.spikeTimes(idx));
            tool.GlobalIdx = idx(order);                 % time-ordered globals

            st = s.spikeTimes(tool.GlobalIdx);
            isi = diff(st) / s.samplingRate * 1000;      % ms, length n-1
            tool.IsiNextMs = [isi; NaN];

            wf = double(s.waveforms(tool.GlobalIdx, :));
            if size(wf, 1) >= 3 && ~isempty(wf)
                [~, sc] = pca(wf, NumComponents=min(3, size(wf, 2)));
                if size(sc, 2) < 3
                    sc(:, end+1:3) = 0;
                end
            else
                sc = zeros(numel(tool.GlobalIdx), 3);
            end
            tool.Scores = sc;

            tool.IsiMask = false(numel(tool.GlobalIdx), 1);
            tool.PCMask = false(numel(tool.GlobalIdx), 1);
        end

        function buildUI(tool)
            tool.UIFigure = uifigure(Name=sprintf("ISI split - cluster %d", ...
                tool.ClusterId), Position=[160 160 1000 560]);
            g = uigridlayout(tool.UIFigure, [2 4]);
            g.RowHeight = {"1x", 40};
            g.ColumnWidth = {"1x", "1x", "1x", "1x"};

            tool.ISIAxes = uiaxes(g);
            tool.ISIAxes.Layout.Row = 1; tool.ISIAxes.Layout.Column = [1 2];

            tool.PCAxes = uiaxes(g);
            tool.PCAxes.Layout.Row = 1; tool.PCAxes.Layout.Column = [3 4];

            circleBtn = uibutton(g, Text="Circle ISI bars", ...
                ButtonPushedFcn=@(~, ~) tool.circleISI());
            circleBtn.Layout.Row = 2; circleBtn.Layout.Column = 1;

            lassoBtn = uibutton(g, Text="Lasso PC features", ...
                ButtonPushedFcn=@(~, ~) tool.lassoPC());
            lassoBtn.Layout.Row = 2; lassoBtn.Layout.Column = 2;

            tool.PCDropdown = uidropdown(g, ...
                Items=["PC1 vs PC2", "PC1 vs PC3", "PC2 vs PC3"], ...
                ValueChangedFcn=@(~, ~) tool.refreshPC());
            tool.PCDropdown.Layout.Row = 2; tool.PCDropdown.Layout.Column = 3;

            splitBtn = uibutton(g, Text="Split off selected", ...
                BackgroundColor=[1 0.85 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.splitSelected());
            splitBtn.Layout.Row = 2; splitBtn.Layout.Column = 4;

            tool.InfoLabel = uilabel(tool.UIFigure, Text="", ...
                Position=[10 4 980 16], FontColor=[0.2 0.2 0.2]);
            tool.setInfo("Circle the ISI bars near 0 to highlight suspect spikes.");
        end

        function circleISI(tool)
            %CIRCLEISI Freehand-draw over ISI bars; select that ISI range.
            roi = drawfreehand(tool.ISIAxes, Color=[0.85 0.5 0.1]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 3
                return
            end
            xRange = [min(roi.Position(:, 1)), max(roi.Position(:, 1))];
            delete(roi);
            tool.selectIsiRange(xRange);
        end

        function lassoPC(tool)
            %LASSOPC Freehand-draw a lasso in the PC scatter to select spikes.
            if ~any(tool.IsiMask)
                tool.setInfo("Circle ISI bars first, then lasso to refine.");
            end
            roi = drawfreehand(tool.PCAxes, Color=[0.8 0.1 0.1]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 3
                return
            end
            poly = roi.Position;
            delete(roi);
            tool.selectPCLasso(poly);
        end

        function refreshISI(tool)
            ax = tool.ISIAxes;
            cla(ax);
            hold(ax, "on");
            edges = 0:tool.BinMs:tool.MaxMs;
            isi = tool.IsiNextMs(~isnan(tool.IsiNextMs));
            histogram(ax, isi, edges, FaceColor=[0.3 0.5 0.8], EdgeColor="none");
            % Overlay the ISIs that belong to the currently highlighted spikes.
            selIsi = tool.IsiNextMs(tool.IsiMask & ~isnan(tool.IsiNextMs));
            if ~isempty(selIsi)
                histogram(ax, selIsi, edges, FaceColor=[0.85 0.5 0.1], ...
                    EdgeColor="none");
            end
            xline(ax, tool.RefractoryMs, "--", Color=[0.6 0 0]);
            hold(ax, "off");
            pct = 100 * sum(isi < tool.RefractoryMs) / max(1, numel(isi));
            title(ax, sprintf("ISI histogram  (%.2f%% < %g ms)", pct, tool.RefractoryMs));
            xlabel(ax, "ISI (ms)");
            ylabel(ax, "Count");
            xlim(ax, [0 tool.MaxMs]);
        end

        function refreshPC(tool)
            ax = tool.PCAxes;
            cla(ax);
            hold(ax, "on");
            [ax1, ax2] = tool.pcColumns();
            x = tool.Scores(:, ax1);
            y = tool.Scores(:, ax2);
            scatter(ax, x, y, 6, [0.6 0.6 0.6], "filled", MarkerFaceAlpha=0.3);
            if any(tool.IsiMask)
                scatter(ax, x(tool.IsiMask), y(tool.IsiMask), 10, ...
                    [0.85 0.5 0.1], "filled", MarkerFaceAlpha=0.6);
            end
            if any(tool.PCMask)
                scatter(ax, x(tool.PCMask), y(tool.PCMask), 14, ...
                    [0.8 0.1 0.1], "filled");
            end
            hold(ax, "off");
            labels = tool.PCDropdown.Value;
            title(ax, "PC features  (" + labels + ")");
            parts = split(labels, " vs ");
            xlabel(ax, parts(1));
            ylabel(ax, parts(2));
            axis(ax, "tight");
        end

        function [ax1, ax2] = pcColumns(tool)
            switch tool.PCDropdown.Value
                case "PC1 vs PC3"
                    ax1 = 1; ax2 = 3;
                case "PC2 vs PC3"
                    ax1 = 2; ax2 = 3;
                otherwise            % "PC1 vs PC2"
                    ax1 = 1; ax2 = 2;
            end
        end

        function setInfo(tool, msg)
            if ~isempty(tool.InfoLabel) && isvalid(tool.InfoLabel)
                tool.InfoLabel.Text = string(msg);
            end
        end
    end
end
