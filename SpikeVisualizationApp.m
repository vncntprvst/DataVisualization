classdef SpikeVisualizationApp < handle
    %SPIKEVISUALIZATIONAPP Interactive viewer and curation tool for sorted spikes.
    %   SpikeVisualizationApp opens a file dialog to choose a sorting to load.
    %
    %   SpikeVisualizationApp(phyDir) opens a Phy/KiloSort export folder
    %   (the folder containing params.py), for example:
    %       SpikeVisualizationApp("D:\...\S113L4A5_12540\tdc\phy_postmerge")
    %
    %   SpikeVisualizationApp(spikes) uses an already-loaded model struct
    %   (see LoadSpikesPhy).
    %
    %   For the selected cluster(s) the window shows overlaid and mean waveforms,
    %   amplitude drift, the (log) inter-spike-interval histogram, the
    %   autocorrelogram, a PCA feature scatter, and a scrollable raw-data excerpt
    %   with spike markers. Voltages are shown in microvolts when a gain is known.
    %
    %   Raw-data excerpt:
    %     * Ctrl + mouse wheel zooms the time axis; the wheel alone pans it.
    %       Left/Right arrows pan, Up/Down arrows zoom.
    %     * "Add threshold" drops a draggable horizontal red line. Every spike
    %       whose waveform crosses all threshold lines is drawn red on the
    %       waveform and mean panels and red-outlined in the amplitude and PC
    %       feature panels. Delete a line with the Delete key or "Clear thresholds".
    %
    %   Clusters can be labelled, merged, saved, and split/curated via the
    %   "Curation" button (see CurationTool).
    %
    %   This is the modern uifigure replacement for the legacy GUIDE
    %   SpikeVisualizationGUI.
    %
    %   See also LoadSpikesPhy, CurationTool.

    properties (SetAccess = private)
        Spikes struct = struct()   % normalized data model (see LoadSpikesPhy)
        Classification table       % one row per labelled cluster
    end

    properties (Access = private)
        UIFigure               matlab.ui.Figure
        ClusterTable           matlab.ui.control.Table
        ClusterRowIds double = []   % cluster id per row of ClusterTable
        ClusterSortBy string = "id" % "id" | "count" | "label"
        ClassTable             matlab.ui.control.Table
        WaveformAxes           matlab.ui.control.UIAxes
        MeanAxes               matlab.ui.control.UIAxes
        AmplitudeAxes          matlab.ui.control.UIAxes
        ISIAxes                matlab.ui.control.UIAxes
        FeatureAxes            matlab.ui.control.UIAxes
        TraceAxes              matlab.ui.control.UIAxes
        TraceSlider            matlab.ui.control.Slider
        InfoLabel              matlab.ui.control.Label
        CurateTool                         % CurationTool handle (companion)

        Traces = []            % memmapfile of the raw data, lazily opened
        TraceCenter double = NaN
        TraceHalfWidth double = NaN
        TimeWindow double = []  % [start stop] s that curation is restricted to

        ThresholdY double = []             % threshold line levels, microvolts
        ThresholdHandles = images.roi.Line.empty
        SuppressLineEvents logical = false
        HighlightIdx logical = logical.empty(0, 1)   % spikes crossing all lines
        ThresholdScope string = "window"  % "window" or "all"
        ScopeMenuWindow matlab.ui.container.Menu
        ScopeMenuAll matlab.ui.container.Menu
        TraceYLim double = []              % fixed voltage limits for the trace
        UpdatingTrace logical = false      % guards the trace XLim listener
        RealignMode string = "peak"        % "peak" | "trough" | "firstpeak"
        RealignMenus struct = struct()     % checkable realign-target menu items

        CorrLayout             matlab.ui.container.GridLayout   % CCG matrix host
        UndoMenu               matlab.ui.container.Menu
        RedoMenu               matlab.ui.container.Menu
        UndoStack cell = {}
        RedoStack cell = {}
    end

    properties (Constant, Access = private)
        MaxWaveformsToPlot = 200
        ISIMinMs = 0.2         % log-axis floor for the ISI histogram
        ISIHistMaxMs = 1000    % log-axis ceiling for the ISI histogram
        ISIMaxMs = 50          % autocorrelogram half-span (ms)
        ISIBinMs = 0.5
        RefractoryMs = 2
        HighlightColor = [0.85 0 0]
        MaxCorrClusters = 4    % cap on the CCG matrix size
        PrefGroup = "SpikeVisualizationApp"
    end

    methods
        function app = SpikeVisualizationApp(source)
            arguments
                source = []    % [] -> ask; folder string -> Phy; struct -> model
            end
            app.Classification = emptyClassificationTable();
            addonsDir = fullfile(fileparts(mfilename("fullpath")), "addons");
            if isfolder(addonsDir)
                addpath(addonsDir);
            end
            app.buildUI();
            try
                if isempty(source)
                    app.promptForSource();
                elseif isstruct(source)
                    app.setData(source);
                else
                    app.loadPhyFolder(string(source));
                end
            catch err
                app.setInfo("Load failed: " + err.message);
            end
        end

        function delete(app)
            if ~isempty(app.CurateTool) && isvalid(app.CurateTool)
                delete(app.CurateTool);
            end
            if isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end

        function f = figureHandle(app)
            %FIGUREHANDLE Return the underlying uifigure (for scripting/testing).
            f = app.UIFigure;
        end
    end

    % ---------------------------------------------------------------- loading
    methods
        function loadPhyFolder(app, phyDir)
            %LOADPHYFOLDER Load a Phy export folder and display its first cluster.
            arguments
                app
                phyDir (1,1) string {mustBeFolder}
            end
            app.setInfo("Loading " + phyDir + " ...");
            drawnow;
            model = LoadSpikesPhy(phyDir);
            app.setData(model);
            setpref(app.PrefGroup, "lastDir", char(phyDir));
        end

        function setData(app, model)
            %SETDATA Adopt a normalized spike model and refresh every view.
            arguments
                app
                model (1,1) struct
            end
            app.Spikes = normalizeModel(model);
            app.Traces = [];
            app.TraceCenter = NaN;
            app.TraceYLim = [];
            app.TimeWindow = [];
            app.ThresholdY = [];
            app.UndoStack = {};
            app.RedoStack = {};
            app.HighlightIdx = false(numel(app.Spikes.clusters), 1);
            app.Classification = emptyClassificationTable();
            app.seedClassificationFromModel();
            app.syncClassification();
            app.loadSavedLabels();      % restore SU/MU/Noise/notes from a prior save
            app.refreshClassTable();
            app.updateClusterList();
            ids = app.Spikes.clusterIds;
            counts = arrayfun(@(c) sum(app.Spikes.clusters == c), ids);
            [~, best] = max(counts);
            app.setSelectedClusters(ids(best));
            app.refreshAll();
            app.setInfo(sprintf("%s  |  %d clusters, %d spikes, %.0f kHz, %.0f s", ...
                app.Spikes.name, numel(ids), numel(app.Spikes.clusters), ...
                app.Spikes.samplingRate / 1000, ...
                app.Spikes.numSamples / app.Spikes.samplingRate));
        end
    end

    % ------------------------------------------------------------------ views
    methods (Access = private)
        function refreshAll(app)
            app.plotTrace();            % sets the trace window, needed for scope
            app.recomputeHighlight();
            app.plotWaveforms();
            app.plotMean();
            app.plotAmplitude();
            app.plotISI();
            app.plotCorrelograms();
            app.plotFeatures();
        end

        function refreshHighlightPanels(app)
            % Re-plot only the panels that show the threshold highlight, leaving
            % the raw-trace axes (and its ROI lines) untouched.
            app.plotWaveforms();
            app.plotMean();
            app.plotAmplitude();
            app.plotFeatures();
        end

        function cids = selectedClusters(app)
            rows = app.ClusterTable.Selection;
            if isempty(rows) || isempty(app.ClusterRowIds)
                cids = double.empty(0, 1);
            else
                cids = double(app.ClusterRowIds(rows(:)));
            end
        end

        function setSelectedClusters(app, ids)
            [tf, rows] = ismember(ids, app.ClusterRowIds);
            app.ClusterTable.Selection = reshape(rows(tf), 1, []);
        end

        function s = voltScale(app)
            % Microvolts per ADC count (1 when no usable gain is known).
            s = 1;
            if isfield(app.Spikes, "uvPerADC")
                s = app.Spikes.uvPerADC;
            end
        end

        function u = voltUnit(app)
            if app.voltScale() == 1
                u = "ADC";
            else
                u = "\muV";
            end
        end

        function hi = highlightSubset(app, globalMask)
            % Indices (subsampled) of highlighted spikes within globalMask.
            hi = find(globalMask & app.HighlightIdx);
            if numel(hi) > app.MaxWaveformsToPlot
                hi = hi(round(linspace(1, numel(hi), app.MaxWaveformsToPlot)));
            end
        end

        function plotWaveforms(app)
            ax = app.WaveformAxes;
            cla(ax);
            hold(ax, "on");
            wf = app.Spikes.waveforms;
            t = app.waveformTimeMs();
            sc = app.voltScale();
            nShown = 0;
            nTotal = 0;
            for cid = app.selectedClusters()'
                idx = find(app.Spikes.clusters == cid);
                if isempty(idx) || isempty(wf)
                    continue
                end
                nTotal = nTotal + numel(idx);
                idx = idx(round(linspace(1, numel(idx), ...
                    min(numel(idx), app.MaxWaveformsToPlot))));
                nShown = nShown + numel(idx);
                plot(ax, t, wf(idx, :)' * sc, Color=[app.clusterColor(cid) 0.15]);
            end
            sel = ismember(app.Spikes.clusters, app.selectedClusters());
            hi = app.highlightSubset(sel);
            if ~isempty(hi)
                plot(ax, t, wf(hi, :)' * sc, Color=app.HighlightColor, LineWidth=0.75);
            end
            hold(ax, "off");
            title(ax, sprintf("Waveforms (%d of %d shown)", nShown, nTotal));
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Voltage (" + app.voltUnit() + ")");
            axis(ax, "tight");
        end

        function plotMean(app)
            ax = app.MeanAxes;
            cla(ax);
            hold(ax, "on");
            wf = app.Spikes.waveforms;
            t = app.waveformTimeMs();
            sc = app.voltScale();
            for cid = app.selectedClusters()'
                idx = app.Spikes.clusters == cid;
                if ~any(idx) || isempty(wf)
                    continue
                end
                m = mean(wf(idx, :), 1) * sc;
                sd = std(double(wf(idx, :)), 0, 1) * sc;   % spread, visible band
                c = app.clusterColor(cid);
                fill(ax, [t fliplr(t)], [m - sd fliplr(m + sd)], c, ...
                    EdgeColor="none", FaceAlpha=0.25);
                plot(ax, t, m, Color=c, LineWidth=1.5);
            end
            sel = ismember(app.Spikes.clusters, app.selectedClusters());
            hi = app.highlightSubset(sel);
            if ~isempty(hi)
                plot(ax, t, wf(hi, :)' * sc, Color=[app.HighlightColor 0.5], ...
                    LineWidth=0.5);
            end
            hold(ax, "off");
            title(ax, "Mean \pm SD");
            xlabel(ax, "Time (ms)");
            ylabel(ax, "Voltage (" + app.voltUnit() + ")");
            axis(ax, "tight");
        end

        function plotAmplitude(app)
            ax = app.AmplitudeAxes;
            cla(ax);
            hold(ax, "on");
            fs = app.Spikes.samplingRate;
            sc = app.voltScale();
            for cid = app.selectedClusters()'
                idx = app.Spikes.clusters == cid;
                if ~any(idx)
                    continue
                end
                scatter(ax, app.Spikes.spikeTimes(idx) / fs, ...
                    app.Spikes.amplitudePP(idx) * sc, 4, app.clusterColor(cid), ...
                    "filled", MarkerFaceAlpha=0.3);
            end
            sel = ismember(app.Spikes.clusters, app.selectedClusters());
            hiMask = sel & app.HighlightIdx;
            if any(hiMask)
                scatter(ax, app.Spikes.spikeTimes(hiMask) / fs, ...
                    app.Spikes.amplitudePP(hiMask) * sc, 16, app.HighlightColor, ...
                    "o", LineWidth=0.75);
            end
            hold(ax, "off");
            title(ax, "Amplitude vs time (drift)");
            xlabel(ax, "Time (s)");
            ylabel(ax, "Peak-to-peak (" + app.voltUnit() + ")");
            axis(ax, "tight");
        end

        function plotISI(app)
            ax = app.ISIAxes;
            cla(ax);
            cids = app.selectedClusters();
            if isempty(cids)
                return
            end
            cid = app.dominantCluster(cids);
            isiMs = app.isiMsOf(cid);
            edges = logspace(log10(app.ISIMinMs), log10(app.ISIHistMaxMs), 64);
            histogram(ax, isiMs, edges, FaceColor=app.clusterColor(cid), ...
                EdgeColor="none");
            set(ax, XScale="log");
            xline(ax, app.RefractoryMs, "--", Color=[0.6 0 0]);
            pctViol = 100 * sum(isiMs < app.RefractoryMs) / max(1, numel(isiMs));
            title(ax, sprintf("ISI  (cluster %d)  %.2f%% < %g ms", ...
                cid, pctViol, app.RefractoryMs));
            xlabel(ax, "ISI (ms, log)");
            ylabel(ax, "Count");
            xlim(ax, [app.ISIMinMs app.ISIHistMaxMs]);
        end

        function plotCorrelograms(app)
            %PLOTCORRELOGRAMS ACG (diagonal) / CCG (off-diagonal) matrix for the
            %   selected clusters. Row i is the reference (lag = t_col - t_row),
            %   so reading the transposed cell flips the reference cluster.
            delete(app.CorrLayout.Children);
            cids = app.selectedClusters();
            cids = cids(1:min(numel(cids), app.MaxCorrClusters));
            n = numel(cids);
            if n == 0
                app.CorrLayout.RowHeight = {"1x"};
                app.CorrLayout.ColumnWidth = {"1x"};
                return
            end
            app.CorrLayout.RowHeight = repmat({"1x"}, 1, n);
            app.CorrLayout.ColumnWidth = repmat({"1x"}, 1, n);
            compact = n > 1;
            for i = 1:n
                for j = 1:n
                    ax = uiaxes(app.CorrLayout);
                    ax.Layout.Row = i; ax.Layout.Column = j;
                    [counts, lags] = app.correlogram(cids(i), cids(j));
                    if i == j
                        faceColor = app.clusterColor(cids(i));
                    else
                        faceColor = [0.45 0.45 0.45];
                    end
                    bar(ax, lags, counts, 1, FaceColor=faceColor, EdgeColor="none");
                    xlim(ax, [-app.ISIMaxMs app.ISIMaxMs]);
                    if compact
                        ax.XTick = []; ax.YTick = [];
                        if i == 1
                            title(ax, "c" + cids(j), FontSize=9);
                        end
                        if j == 1
                            ylabel(ax, "c" + cids(i), FontSize=9);
                        end
                    else
                        title(ax, sprintf("ACG (cluster %d)", cids(i)));
                        xlabel(ax, "Lag (ms)");
                        ylabel(ax, "Count");
                    end
                end
            end
        end

        function plotFeatures(app)
            ax = app.FeatureAxes;
            cla(ax);
            hold(ax, "on");
            for cid = app.selectedClusters()'
                idxCid = find(app.Spikes.clusters == cid);
                sc = app.featureScores(cid);
                if isempty(sc)
                    continue
                end
                scatter(ax, sc(:, 1), sc(:, 2), 5, app.clusterColor(cid), ...
                    "filled", MarkerFaceAlpha=0.3);
                hiWithin = app.HighlightIdx(idxCid);
                if any(hiWithin)
                    scatter(ax, sc(hiWithin, 1), sc(hiWithin, 2), 16, ...
                        app.HighlightColor, "o", LineWidth=0.75);
                end
            end
            hold(ax, "off");
            title(ax, "PC features (per-cluster PCA of waveforms)");
            xlabel(ax, "PC1");
            ylabel(ax, "PC2");
            axis(ax, "tight");
        end

        function plotTrace(app)
            app.UpdatingTrace = true;   % suppress the XLim listener while we draw
            ax = app.TraceAxes;
            cla(ax);
            if ~app.ensureTraces()
                title(ax, "Raw trace unavailable (no data file)");
                app.UpdatingTrace = false;
                return
            end
            fs = app.Spikes.samplingRate;
            sc = app.voltScale();
            app.ensureTraceCenter();
            [first, last] = app.traceWindowBounds();
            ch = app.Spikes.dataChannel + 1;
            excerpt = single(app.Traces.Data.raw(ch, first:last)) * sc;
            tSec = (first:last) / fs;
            plot(ax, tSec, excerpt, Color=[0 0 0 0.8], LineWidth=0.1);
            hold(ax, "on");
            yMarker = app.TraceYLim(1) + 0.03 * diff(app.TraceYLim);
            for cid = app.selectedClusters()'
                st = app.Spikes.spikeTimes(app.Spikes.clusters == cid);
                st = st(st >= first & st <= last);
                plot(ax, st / fs, repmat(yMarker, size(st)), "^", ...
                    Color=app.clusterColor(cid), MarkerFaceColor=app.clusterColor(cid), ...
                    MarkerSize=4, LineStyle="none");
            end
            hold(ax, "off");
            title(ax, sprintf("Raw data excerpt  (%.0f ms window)", ...
                2 * app.TraceHalfWidth / fs * 1000));
            xlabel(ax, "Time (s)");
            ylabel(ax, "Voltage (" + app.voltUnit() + ")");
            xlim(ax, [tSec(1) tSec(end)]);
            ylim(ax, app.TraceYLim);        % fixed across pan/zoom
            app.TraceSlider.Value = min(1, max(0, ...
                app.TraceCenter / app.Spikes.numSamples));   % keep slider in sync
            app.drawThresholdLines();
            app.UpdatingTrace = false;
        end

        function onTraceXLimChanged(app)
            %ONTRACEXLIMCHANGED Reload the excerpt and respan threshold lines when
            %   the trace x-limits change by any means (e.g. the axes zoom/pan
            %   toolbar), not only via the app's own zoom.
            if app.UpdatingTrace || isnan(app.TraceCenter) || ~app.ensureTraces()
                return
            end
            fs = app.Spikes.samplingRate;
            xl = xlim(app.TraceAxes);
            first = max(1, round(xl(1) * fs));
            last = min(app.Spikes.numSamples, round(xl(2) * fs));
            if last - first < 2
                return
            end
            app.TraceCenter = round((first + last) / 2);
            app.TraceHalfWidth = max(1, round((last - first) / 2));
            app.updateTrace();
        end

        function ensureTraceCenter(app)
            if isnan(app.TraceCenter)
                fs = app.Spikes.samplingRate;
                app.TraceHalfWidth = round(fs / 2);            % 1 s window
                app.TraceCenter = round(app.Spikes.numSamples / 2);   % mid-recording
            end
        end

        function [first, last] = traceWindowBounds(app)
            app.ensureTraceCenter();
            half = app.TraceHalfWidth;
            first = max(1, round(app.TraceCenter) - half);
            last = min(app.Spikes.numSamples, round(app.TraceCenter) + half);
        end

        function updateTrace(app)
            % Redraw the trace and, when thresholds are scoped to the visible
            % window, recompute which spikes are highlighted.
            app.plotTrace();
            if ~isempty(app.ThresholdY) && app.ThresholdScope == "window"
                app.recomputeHighlight();
                app.refreshHighlightPanels();
            end
        end
    end

    % ------------------------------------------------------------ computations
    methods (Access = private)
        function t = waveformTimeMs(app)
            w = app.Spikes.waveformWindow;
            t = (-w(1):w(2)) / app.Spikes.samplingRate * 1000;
            if numel(t) ~= size(app.Spikes.waveforms, 2)
                t = 1:size(app.Spikes.waveforms, 2);
            end
        end

        function isiMs = isiMsOf(app, cid)
            st = sort(app.Spikes.spikeTimes(app.Spikes.clusters == cid));
            isiMs = diff(st) / app.Spikes.samplingRate * 1000;
        end

        function [counts, lags] = correlogram(app, cidA, cidB)
            fs = app.Spikes.samplingRate;
            a = sort(app.Spikes.spikeTimes(app.Spikes.clusters == cidA)) / fs * 1000;
            b = sort(app.Spikes.spikeTimes(app.Spikes.clusters == cidB)) / fs * 1000;
            % Cap the reference train for responsiveness (shape is preserved).
            maxRef = 15000;
            if numel(a) > maxRef
                a = a(round(linspace(1, numel(a), maxRef)));
            end
            edges = -app.ISIMaxMs:app.ISIBinMs:app.ISIMaxMs;
            lags = edges(1:end-1) + app.ISIBinMs / 2;
            counts = zeros(1, numel(lags));
            for k = 1:numel(a)
                d = b - a(k);
                d = d(abs(d) <= app.ISIMaxMs);
                counts = counts + histcounts(d, edges);
            end
            if cidA == cidB
                counts(abs(lags) < app.ISIBinMs) = 0;   % drop the zero-lag peak
            end
        end

        function sc = featureScores(app, cid)
            %FEATURESCORES PCA scores (first components) of one cluster's waveforms.
            idx = app.Spikes.clusters == cid;
            wf = double(app.Spikes.waveforms(idx, :));
            if size(wf, 1) < 3 || isempty(wf)
                sc = [];
                return
            end
            [~, sc] = pca(wf, NumComponents=min(3, size(wf, 2)));
        end

        function cid = dominantCluster(app, cids)
            counts = arrayfun(@(c) sum(app.Spikes.clusters == c), cids);
            [~, k] = max(counts);
            cid = cids(k);
        end

        function pos = dominantFirstSpike(app)
            cid = app.dominantCluster(app.selectedClusters());
            pos = find(app.Spikes.clusters == cid, 1);
            if isempty(pos)
                pos = 1;
            end
        end

        function ok = ensureTraces(app)
            ok = false;
            if ~isempty(app.Traces)
                ok = true;
                return
            end
            if ~isfield(app.Spikes, "datPath") || ~isfile(app.Spikes.datPath) ...
                    || isnan(app.Spikes.numSamples)
                return
            end
            app.Traces = memmapfile(app.Spikes.datPath, ...
                Format={'int16', [app.Spikes.numChannels app.Spikes.numSamples], 'raw'});
            app.computeTraceYLim();
            ok = true;
        end

        function computeTraceYLim(app)
            % Robust, fixed voltage limits sampled across the whole recording so
            % the trace axis does not rescale while panning.
            sc = app.voltScale();
            ch = app.Spikes.dataChannel + 1;
            stride = max(1, floor(app.Spikes.numSamples / 500000));
            sample = single(app.Traces.Data.raw(ch, 1:stride:end)) * sc;
            lim = prctile(double(sample), [0.05 99.95]);
            pad = 0.15 * diff(lim);
            app.TraceYLim = [lim(1) - pad, lim(2) + pad];
        end

    end

    % ------------------------------------------------------------ trace zoom/pan
    methods (Access = private)
        function onScroll(app, ~, event)
            if ~app.ensureTraces() || isnan(app.TraceCenter)
                return
            end
            if any(strcmp(app.UIFigure.CurrentModifier, "control"))
                app.zoomTrace(1.25 ^ event.VerticalScrollCount);
            else
                app.panTrace(0.2 * event.VerticalScrollCount);
            end
        end

        function onKey(app, ~, event)
            if ~app.ensureTraces() || isnan(app.TraceCenter)
                return
            end
            switch event.Key
                case "leftarrow"
                    app.panTrace(-0.25);
                case "rightarrow"
                    app.panTrace(0.25);
                case {"uparrow", "add", "equal"}
                    app.zoomTrace(1 / 1.25);
                case {"downarrow", "subtract", "hyphen"}
                    app.zoomTrace(1.25);
                otherwise
                    % other keys are not handled by the trace view
            end
        end

        function zoomTrace(app, factor)
            fs = app.Spikes.samplingRate;
            minHalf = round(fs * 0.005);
            maxHalf = round(min(app.Spikes.numSamples / 2, fs * 10));
            app.TraceHalfWidth = min(maxHalf, max(minHalf, ...
                round(app.TraceHalfWidth * factor)));
            app.clampTraceCenter();
            app.updateTrace();
        end

        function panTrace(app, fractionOfWindow)
            app.TraceCenter = app.TraceCenter + ...
                fractionOfWindow * 2 * app.TraceHalfWidth;
            app.clampTraceCenter();
            app.TraceSlider.Value = min(1, max(0, ...
                app.TraceCenter / app.Spikes.numSamples));
            app.updateTrace();
        end

        function clampTraceCenter(app)
            half = app.TraceHalfWidth;
            app.TraceCenter = min(app.Spikes.numSamples - half, ...
                max(half + 1, app.TraceCenter));
        end
    end

    % --------------------------------------------------- threshold-line highlight
    methods (Access = private)
        function addThresholdLine(app)
            if ~app.ensureTraces()
                app.setInfo("Load a recording with raw data to add thresholds.");
                return
            end
            yl = ylim(app.TraceAxes);
            app.ThresholdY(end + 1) = mean(yl);
            app.drawThresholdLines();
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end

        function clearThresholdLines(app)
            app.ThresholdY = [];
            app.drawThresholdLines();
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end

        function drawThresholdLines(app)
            app.SuppressLineEvents = true;
            delete(app.ThresholdHandles(isvalid(app.ThresholdHandles)));
            app.ThresholdHandles = images.roi.Line.empty;
            if ~isempty(app.ThresholdY) && isvalid(app.TraceAxes)
                xl = xlim(app.TraceAxes);
                for k = 1:numel(app.ThresholdY)
                    y = app.ThresholdY(k);
                    roi = images.roi.Line(app.TraceAxes, ...
                        Position=[xl(1) y; xl(2) y], Color=[0.85 0 0], ...
                        LineWidth=1, MarkerSize=0.1);
                    addlistener(roi, "MovingROI", @(s, ~) app.constrainLine(s));
                    addlistener(roi, "ROIMoved", @(~, ~) app.onLineMoved());
                    addlistener(roi, "DeletingROI", @(s, ~) app.onLineDeleted(s));
                    app.ThresholdHandles(k) = roi;
                end
            end
            app.SuppressLineEvents = false;
        end

        function constrainLine(app, roi)
            if app.SuppressLineEvents
                return
            end
            ym = mean(roi.Position(:, 2));
            xl = xlim(app.TraceAxes);
            roi.Position = [xl(1) ym; xl(2) ym];   % keep the line horizontal
        end

        function onLineMoved(app)
            if app.SuppressLineEvents
                return
            end
            app.syncThresholdY();
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end

        function onLineDeleted(app, roi)
            if app.SuppressLineEvents
                return
            end
            keep = arrayfun(@(h) isvalid(h) && h ~= roi, app.ThresholdHandles);
            app.ThresholdHandles = app.ThresholdHandles(keep);
            app.syncThresholdY();
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end

        function syncThresholdY(app)
            valid = app.ThresholdHandles(isvalid(app.ThresholdHandles));
            app.ThresholdY = arrayfun(@(h) mean(h.Position(:, 2)), valid);
        end

        function recomputeHighlight(app)
            app.HighlightIdx = false(numel(app.Spikes.clusters), 1);
            if isempty(app.ThresholdY) || ~isfield(app.Spikes, "wfMinADC") ...
                    || isempty(app.Spikes.wfMinADC)
                return
            end
            sc = app.voltScale();
            loADC = min(app.ThresholdY) / sc;
            hiADC = max(app.ThresholdY) / sc;
            sel = ismember(app.Spikes.clusters, app.selectedClusters());
            crossesAll = app.Spikes.wfMinADC < loADC & app.Spikes.wfMaxADC > hiADC;
            mask = sel & crossesAll;
            if app.ThresholdScope == "window" && ~isnan(app.TraceCenter)
                [first, last] = app.traceWindowBounds();
                mask = mask & app.Spikes.spikeTimes >= first ...
                    & app.Spikes.spikeTimes <= last;
            end
            app.HighlightIdx = mask;
        end
    end

    methods
        function setThresholds(app, yMicrovolts)
            %SETTHRESHOLDS Set threshold levels programmatically (scriptable/testing).
            %   Equivalent to adding horizontal threshold lines at the given
            %   microvolt levels on the raw-trace panel.
            arguments
                app
                yMicrovolts (1, :) double
            end
            app.ThresholdY = yMicrovolts;
            app.drawThresholdLines();
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end
    end

    % ---------------------------------------------------------- classification
    methods (Access = private)
        function classifySelected(app, label)
            app.pushUndo();
            cids = app.selectedClusters();
            for cid = cids'
                row = app.Classification.ClusterID == cid;
                if any(row)
                    app.Classification.Label(row) = string(label);
                else
                    app.Classification = [app.Classification; ...
                        {cid, string(label), ""}];
                end
            end
            app.Classification = sortrows(app.Classification, "ClusterID");
            app.refreshClassTable();
            app.updateClusterList();
        end

        function seedClassificationFromModel(app)
            if ~isfield(app.Spikes, "clusterLabels")
                return
            end
            lbl = app.Spikes.clusterLabels;
            for k = 1:numel(app.Spikes.clusterIds)
                app.Classification = [app.Classification; ...
                    {app.Spikes.clusterIds(k), lbl(k), ""}];
            end
        end

        function refreshClassTable(app)
            app.ClassTable.Data = app.Classification;
        end

        function loadSavedLabels(app)
            %LOADSAVEDLABELS Merge labels/notes from a prior cluster_curation.csv.
            if ~isfield(app.Spikes, "phyDir")
                return
            end
            csv = fullfile(app.Spikes.phyDir, "cluster_curation.csv");
            if ~isfile(csv)
                return
            end
            try
                saved = readtable(csv, TextType="string");
            catch
                return   % unreadable/legacy file: keep the Phy group labels
            end
            vars = string(saved.Properties.VariableNames);
            if ~ismember("ClusterID", vars)
                return
            end
            for i = 1:height(saved)
                row = app.Classification.ClusterID == saved.ClusterID(i);
                if ~any(row)
                    continue
                end
                if ismember("Label", vars)
                    app.Classification.Label(row) = missingToEmpty(saved.Label(i));
                end
                if ismember("Note", vars)
                    app.Classification.Note(row) = missingToEmpty(saved.Note(i));
                end
            end
        end

        function syncClassification(app)
            %SYNCCLASSIFICATION Keep one label row per existing cluster so that
            %   newly split clusters appear in the table right away.
            ids = app.Spikes.clusterIds;
            keep = ismember(app.Classification.ClusterID, ids);
            app.Classification = app.Classification(keep, :);
            missing = setdiff(ids, app.Classification.ClusterID);
            for cid = missing(:)'
                app.Classification = [app.Classification; {cid, "", ""}];
            end
            app.Classification = sortrows(app.Classification, "ClusterID");
            app.refreshClassTable();
        end

        function onClassEdited(app, event)
            app.pushUndo();
            r = event.Indices(1);
            c = event.Indices(2);
            varName = app.ClassTable.ColumnName{c};
            app.Classification.(varName)(r) = string(event.NewData);
            if varName == "Label"
                app.updateClusterList();
            end
        end

        function saveClassification(app)
            if isempty(app.Spikes) || ~isfield(app.Spikes, "phyDir")
                app.setInfo("Nothing to save.");
                return
            end
            outFile = fullfile(app.Spikes.phyDir, "cluster_curation.csv");
            writetable(app.Classification, outFile);
            app.setInfo("Saved labels to " + outFile);
        end

        function saveSorting(app)
            if ~isfield(app.Spikes, "phyDir") || app.Spikes.sourceType ~= "phy"
                app.setInfo("Save sorting only supported for Phy sources.");
                return
            end
            % All per-spike arrays are written together in one time-sorted order
            % so they stay aligned (e.g. after a realign changes the order).
            [times, order] = sort(app.Spikes.spikeTimes);
            targets = "spike_times.npy";
            datas = {int64(round(times))};
            targets(end + 1) = "spike_clusters.npy";
            datas{end + 1} = int32(app.Spikes.clusters(order));
            if isfield(app.Spikes, "templates") && ~isempty(app.Spikes.templates)
                targets(end + 1) = "spike_templates.npy";
                datas{end + 1} = int32(app.Spikes.templates(order));
            end
            if isfield(app.Spikes, "amplitudes") && ~isempty(app.Spikes.amplitudes)
                targets(end + 1) = "amplitudes.npy";
                datas{end + 1} = double(app.Spikes.amplitudes(order));
            end
            if isfield(app.Spikes, "pcFeatures") && ~isempty(app.Spikes.pcFeatures)
                targets(end + 1) = "pc_features.npy";
                datas{end + 1} = app.Spikes.pcFeatures(order, :, :);
            end
            try
                app.saveArraysAtomic(targets, datas);
            catch err
                app.setInfo("Save failed: " + err.message);
                return
            end
            app.saveClassification();
            app.setInfo(sprintf("Saved sorting (%d arrays) + labels; originals " + ...
                "backed up as *.bak.npy.", numel(targets)));
        end

        function saveArraysAtomic(app, targets, datas)
            %SAVEARRAYSATOMIC Write .npy arrays all-or-nothing (temp then move,
            %   with per-save rollback), so a locked/failed write never leaves the
            %   folder half-updated. Keeps a pristine one-time *.bak.npy too.
            phyDir = app.Spikes.phyDir;
            tmp = fullfile(phyDir, targets + ".savetmp");
            roll = fullfile(phyDir, targets + ".savebak");
            cleaner = onCleanup(@() app.deleteFiles([tmp, roll]));

            for i = 1:numel(targets)
                writeNPY(datas{i}, tmp(i));   % new files; throws if folder unwritable
            end

            moved = false(1, numel(targets));
            try
                for i = 1:numel(targets)
                    tgt = fullfile(phyDir, targets(i));
                    app.backupOriginalOnce(targets(i));
                    if isfile(tgt)
                        copyfile(tgt, roll(i));   % per-save rollback copy
                    end
                    app.clearReadOnly(tgt);
                    app.moveWithRetry(tmp(i), tgt);
                    moved(i) = true;
                end
            catch err
                for i = find(moved)
                    if isfile(roll(i))
                        app.clearReadOnly(fullfile(phyDir, targets(i)));
                        copyfile(roll(i), fullfile(phyDir, targets(i)));
                    end
                end
                error("SpikeVisualizationApp:saveLocked", ...
                    "A sorting file is locked, so nothing was changed (folder " + ...
                    "left intact). Close Phy and any file-sync/antivirus using " + ...
                    "this folder, then retry. Details: %s", err.message);
            end
        end

        function backupOriginalOnce(app, fileName)
            [~, base, ext] = fileparts(fileName);
            src = fullfile(app.Spikes.phyDir, fileName);
            backup = fullfile(app.Spikes.phyDir, base + ".bak" + ext);
            if isfile(src) && ~isfile(backup)
                copyfile(src, backup);
            end
        end

        function clearReadOnly(~, file)
            if isfile(file)
                try
                    fileattrib(file, "+w");
                catch
                    % best effort; a genuine write failure is handled by caller
                end
            end
        end

        function moveWithRetry(app, src, dst)
            app.clearReadOnly(dst);
            lastErr = [];
            for attempt = 1:6
                try
                    movefile(src, dst);
                    return
                catch err
                    lastErr = err;
                    pause(0.4 * attempt);
                    app.clearReadOnly(dst);
                end
            end
            rethrow(lastErr);
        end

        function deleteFiles(~, files)
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
    end

    % ---------------------------------------------- public API for curate tool
    methods
        function cid = currentCluster(app)
            %CURRENTCLUSTER The single most-numerous selected cluster id.
            cids = app.selectedClusters();
            if isempty(cids)
                cid = [];
            else
                cid = app.dominantCluster(cids);
            end
        end

        function cids = selectedClusterIds(app)
            %SELECTEDCLUSTERIDS Cluster ids currently selected in the list.
            cids = app.selectedClusters();
        end

        function [tSec, uv] = traceExcerpt(app, centerSec, halfSec)
            %TRACEEXCERPT Data-channel excerpt (seconds, microvolts) around a time.
            [tSec, uv] = deal([]);
            if ~app.ensureTraces()
                return
            end
            fs = app.Spikes.samplingRate;
            half = round(halfSec * fs);
            c = round(centerSec * fs);
            first = max(1, c - half);
            last = min(app.Spikes.numSamples, c + half);
            ch = app.Spikes.dataChannel + 1;
            raw = single(app.Traces.Data.raw(ch, first:last));
            scv = 1;
            if isfield(app.Spikes, "uvPerADC")
                scv = app.Spikes.uvPerADC;
            end
            uv = raw * scv;
            tSec = (first:last) / fs;
        end

        function w = timeWindowSec(app)
            %TIMEWINDOWSEC Current curation time window [start stop] in seconds.
            dur = app.Spikes.numSamples / app.Spikes.samplingRate;
            if numel(app.TimeWindow) == 2
                w = [max(0, app.TimeWindow(1)), min(dur, app.TimeWindow(2))];
            else
                w = [0, dur];
            end
        end

        function setTimeWindow(app, w)
            %SETTIMEWINDOW Restrict curation to spikes within [start stop] seconds.
            %   Refreshes the open Curation tool. An empty window means "all".
            arguments
                app
                w double
            end
            if numel(w) == 2
                app.TimeWindow = sort(reshape(w, 1, []));
            else
                app.TimeWindow = [];
            end
            if ~isempty(app.CurateTool) && isvalid(app.CurateTool)
                app.CurateTool.refreshWindow();
            end
        end

        function c = clusterColor(app, cid)
            %CLUSTERCOLOR Well-separated display colour for a cluster.
            %   Golden-angle hue spacing over the cluster's position keeps colours
            %   as distinct as possible for any number of clusters (no repeat
            %   after 7 like the default palette).
            k = find(app.Spikes.clusterIds == cid, 1);
            if isempty(k)
                k = 1;
            end
            hue = mod((k - 1) * 0.6180339887, 1);
            c = hsv2rgb([hue, 0.62, 0.92]);
        end

        function applySplit(app, spikeIndices, newCluster)
            %APPLYSPLIT Reassign the given global spike indices to a new cluster id.
            arguments
                app
                spikeIndices (:, 1) double
                newCluster (1, 1) double
            end
            app.pushUndo();
            app.Spikes.clusters(spikeIndices) = newCluster;
            app.Spikes.clusterIds = unique(app.Spikes.clusters);
            app.syncClassification();
            app.updateClusterList();
            app.refreshAll();
            app.setInfo(sprintf("Split %d spikes into cluster %d.", ...
                numel(spikeIndices), newCluster));
        end

        function newId = applyMerge(app, clusterIds)
            %APPLYMERGE Merge the given clusters into their lowest id.
            arguments
                app
                clusterIds (:, 1) double
            end
            newId = [];
            clusterIds = unique(clusterIds);
            if numel(clusterIds) < 2
                app.setInfo("Select at least two clusters to merge.");
                return
            end
            app.pushUndo();
            newId = min(clusterIds);
            app.Spikes.clusters(ismember(app.Spikes.clusters, clusterIds)) = newId;
            app.Spikes.clusterIds = unique(app.Spikes.clusters);
            app.syncClassification();
            app.updateClusterList();
            app.setSelectedClusters(newId);
            app.refreshAll();
            app.setInfo(sprintf("Merged clusters %s into %d.", ...
                num2str(clusterIds'), newId));
        end
    end

    % ------------------------------------------------------------- curate tool
    methods (Access = private)
        function launchCurateTool(app)
            cids = app.selectedClusters();
            if isempty(cids)
                app.setInfo("Select one or more clusters to curate.");
                return
            end
            if ~isempty(app.CurateTool) && isvalid(app.CurateTool)
                delete(app.CurateTool);
            end
            app.CurateTool = CurationTool(app, cids);
        end

        function mergeSelected(app)
            app.applyMerge(app.selectedClusters());
        end

        function realignSelected(app)
            %REALIGNSELECTED Shift each selected cluster's spike times so its mean
            %   waveform peak lands at the sample the spike time is anchored to,
            %   then re-extract the waveforms. Aligns time-shifted duplicate
            %   clusters so they overlay (and can then be merged).
            if ~app.ensureTraces()
                app.setInfo("Realign needs the raw data file.");
                return
            end
            cids = app.selectedClusters();
            if isempty(cids)
                return
            end
            app.pushUndo();
            target = app.Spikes.waveformWindow(1) + 1;   % spike-time sample index
            for cid = cids'
                idx = app.Spikes.clusters == cid;
                if ~any(idx)
                    continue
                end
                m = mean(app.Spikes.waveforms(idx, :), 1);
                shift = app.alignmentSample(m) - target;
                if shift ~= 0
                    app.Spikes.spikeTimes(idx) = app.Spikes.spikeTimes(idx) + shift;
                end
            end
            sel = ismember(app.Spikes.clusters, cids);
            app.Spikes.waveforms(sel, :) = ...
                app.extractWaveformsAt(app.Spikes.spikeTimes(sel));
            app.Spikes.wfMinADC(sel) = double(min(app.Spikes.waveforms(sel, :), [], 2));
            app.Spikes.wfMaxADC(sel) = double(max(app.Spikes.waveforms(sel, :), [], 2));
            app.Spikes.amplitudePP(sel) = app.Spikes.wfMaxADC(sel) - app.Spikes.wfMinADC(sel);
            app.refreshAll();
            app.setInfo(sprintf("Realigned %d cluster(s) to their waveform peak.", ...
                numel(cids)));
        end

        function sampleIdx = alignmentSample(app, meanWaveform)
            %ALIGNMENTSAMPLE Sample of the feature to realign to, per RealignMode.
            switch app.RealignMode
                case "trough"
                    [~, sampleIdx] = min(meanWaveform);
                case "firstpeak"
                    absWf = abs(meanWaveform);
                    [~, locs] = findpeaks(absWf, MinPeakHeight=0.3 * max(absWf));
                    if isempty(locs)
                        [~, sampleIdx] = max(absWf);
                    else
                        sampleIdx = locs(1);
                    end
                otherwise   % "peak" - largest absolute deflection
                    [~, sampleIdx] = max(abs(meanWaveform));
            end
        end

        function wf = extractWaveformsAt(app, times)
            %EXTRACTWAVEFORMSAT Waveforms of the data channel around given samples.
            w = app.Spikes.waveformWindow;
            windowLength = w(1) + w(2) + 1;
            ch = app.Spikes.dataChannel + 1;
            numSamples = app.Spikes.numSamples;
            channelData = single(app.Traces.Data.raw(ch, :));
            times = round(times(:));
            wf = zeros(numel(times), windowLength, "single");
            starts = times - w(1);
            for k = 1:numel(times)
                idx = starts(k):(starts(k) + windowLength - 1);
                valid = idx >= 1 & idx <= numSamples;
                wf(k, valid) = channelData(idx(valid));
            end
        end
    end

    % --------------------------------------------------------------- plumbing
    methods (Access = private)
        function promptForSource(app)
            startDir = getpref(app.PrefGroup, "lastDir", pwd);
            if ~isfolder(startDir)
                startDir = pwd;
            end
            folder = uigetdir(startDir, "Select a Phy export folder (contains params.py)");
            figure(app.UIFigure);
            if isequal(folder, 0)
                app.setInfo("No folder selected. Use File controls to load a sorting.");
                return
            end
            app.loadPhyFolder(string(folder));
        end

        function setClusterSort(app, sortBy)
            app.ClusterSortBy = sortBy;
            app.updateClusterList();
        end

        function updateClusterList(app)
            ids = app.Spikes.clusterIds;
            counts = arrayfun(@(c) sum(app.Spikes.clusters == c), ids);
            labels = strings(numel(ids), 1);
            for k = 1:numel(ids)
                row = app.Classification.ClusterID == ids(k);
                if any(row)
                    labels(k) = app.Classification.Label(find(row, 1));
                end
            end
            switch app.ClusterSortBy   % row order in the table (colours unchanged)
                case "count"
                    [~, ord] = sort(counts, "descend");
                case "label"
                    [~, ord] = sort(labels);
                otherwise
                    ord = 1:numel(ids);
            end
            ids = ids(ord); counts = counts(ord); labels = labels(ord);

            keptSelection = app.selectedClusters();
            data = cell(numel(ids), 4);
            for k = 1:numel(ids)
                data(k, :) = {'', ids(k), counts(k), char(labels(k))};
            end
            app.ClusterRowIds = ids;
            app.ClusterTable.Data = data;
            % Colour swatch in the first column, one per cluster.
            removeStyle(app.ClusterTable);
            for k = 1:numel(ids)
                addStyle(app.ClusterTable, ...
                    uistyle(BackgroundColor=app.clusterColor(ids(k))), ...
                    "cell", [k 1]);
            end
            app.setSelectedClusters(keptSelection);
        end

        function setInfo(app, msg)
            if ~isempty(app.InfoLabel) && isvalid(app.InfoLabel)
                app.InfoLabel.Text = string(msg);
            end
        end

        function onClusterChanged(app, ~, ~)
            app.refreshAll();   % keep the trace where the user left it
        end

        function onTraceSlider(app, ~, event)
            if ~app.ensureTraces()
                return
            end
            app.TraceCenter = event.Value * app.Spikes.numSamples;
            app.clampTraceCenter();
            app.updateTrace();
        end

        function setThresholdScope(app, scope)
            app.ThresholdScope = scope;
            app.ScopeMenuWindow.Checked = matlab.lang.OnOffSwitchState(scope == "window");
            app.ScopeMenuAll.Checked = matlab.lang.OnOffSwitchState(scope == "all");
            app.recomputeHighlight();
            app.refreshHighlightPanels();
        end
    end

    % ----------------------------------------------------------------- buildUI
    methods (Access = private)
        function buildUI(app)
            app.UIFigure = uifigure(Name="SpikeVisualizationApp", ...
                Position=[60 40 1440 1000], ...
                WindowScrollWheelFcn=@(s, e) app.onScroll(s, e), ...
                WindowKeyPressFcn=@(s, e) app.onKey(s, e));
            outer = uigridlayout(app.UIFigure, [2 2]);
            outer.RowHeight = {"1x", 24};
            outer.ColumnWidth = {300, "1x"};

            app.buildMenus();
            app.buildControls(outer);
            app.buildPlots(outer);

            app.InfoLabel = uilabel(outer, Text="Ready.", FontColor=[0.2 0.2 0.2]);
            app.InfoLabel.Layout.Row = 2;
            app.InfoLabel.Layout.Column = [1 2];
        end

        function buildMenus(app)
            edit = uimenu(app.UIFigure, Text="Edit");
            app.UndoMenu = uimenu(edit, Text="Undo", Accelerator="Z", ...
                Enable="off", MenuSelectedFcn=@(~, ~) app.undo());
            app.RedoMenu = uimenu(edit, Text="Redo", Accelerator="Y", ...
                Enable="off", MenuSelectedFcn=@(~, ~) app.redo());

            opt = uimenu(app.UIFigure, Text="Options");
            app.buildRunScriptMenu(opt);
            thr = uimenu(opt, Text="Threshold highlights spikes");
            app.ScopeMenuWindow = uimenu(thr, Text="In current window", ...
                Checked="on", MenuSelectedFcn=@(~, ~) app.setThresholdScope("window"));
            app.ScopeMenuAll = uimenu(thr, Text="In whole recording", ...
                MenuSelectedFcn=@(~, ~) app.setThresholdScope("all"));

            ra = uimenu(opt, Text="Realign to");
            app.RealignMenus.peak = uimenu(ra, Text="Largest peak", Checked="on", ...
                MenuSelectedFcn=@(~, ~) app.setRealignMode("peak"));
            app.RealignMenus.trough = uimenu(ra, Text="Trough (negative)", ...
                MenuSelectedFcn=@(~, ~) app.setRealignMode("trough"));
            app.RealignMenus.firstpeak = uimenu(ra, Text="First peak", ...
                MenuSelectedFcn=@(~, ~) app.setRealignMode("firstpeak"));
        end

        function setRealignMode(app, mode)
            app.RealignMode = mode;
            app.RealignMenus.peak.Checked = ...
                matlab.lang.OnOffSwitchState(mode == "peak");
            app.RealignMenus.trough.Checked = ...
                matlab.lang.OnOffSwitchState(mode == "trough");
            app.RealignMenus.firstpeak.Checked = ...
                matlab.lang.OnOffSwitchState(mode == "firstpeak");
            app.setInfo("Realign target set to: " + mode);
        end

        % --- add-on scripts (Options > Run script) ---
        function buildRunScriptMenu(app, parent)
            run = uimenu(parent, Text="Run script");
            files = app.addonFiles();
            if isempty(files)
                uimenu(run, Text="(no add-on scripts)", Enable="off");
                return
            end
            for k = 1:numel(files)
                f = files(k);
                uimenu(run, Text=addonDisplayName(f), ...
                    MenuSelectedFcn=@(~, ~) app.runAddon(f));
            end
        end

        function files = addonFiles(app) %#ok<MANU>
            dirPath = fullfile(fileparts(mfilename("fullpath")), "addons");
            listing = dir(fullfile(dirPath, "*.m"));
            files = string({listing.folder}') + filesep + string({listing.name}');
            files = files(:);
        end

        function runAddon(app, file)
            if ~isfield(app.Spikes, "phyDir") || isempty(app.Spikes.phyDir)
                app.setInfo("Load a Phy folder before running a script.");
                return
            end
            [~, fn] = fileparts(file);
            app.setInfo("Running " + fn + " ...");
            drawnow;
            app.Traces = [];   % release the raw-data handle so the app never locks the folder
            try
                feval(fn, app.Spikes.phyDir);
                app.reloadData();
                app.setInfo(fn + " finished; data reloaded.");
            catch err
                % On failure the add-on leaves the folder unchanged, so the
                % in-memory data is still valid; the memmap re-acquires lazily.
                app.setInfo(fn + " failed: " + err.message);
            end
        end

        % --- undo / redo ---
        function pushUndo(app)
            if isempty(app.Spikes) || ~isfield(app.Spikes, "clusters")
                return
            end
            app.UndoStack{end + 1} = app.snapshotState();
            if numel(app.UndoStack) > 20
                app.UndoStack(1) = [];
            end
            app.RedoStack = {};
            app.updateUndoMenus();
        end

        function s = snapshotState(app)
            s.clusters = app.Spikes.clusters;
            s.spikeTimes = app.Spikes.spikeTimes;
            s.waveforms = app.Spikes.waveforms;
            s.wfMinADC = app.Spikes.wfMinADC;
            s.wfMaxADC = app.Spikes.wfMaxADC;
            s.amplitudePP = app.Spikes.amplitudePP;
            s.clusterIds = app.Spikes.clusterIds;
            s.classification = app.Classification;
            s.selection = app.selectedClusters();
        end

        function restoreState(app, s)
            app.Spikes.clusters = s.clusters;
            app.Spikes.spikeTimes = s.spikeTimes;
            app.Spikes.waveforms = s.waveforms;
            app.Spikes.wfMinADC = s.wfMinADC;
            app.Spikes.wfMaxADC = s.wfMaxADC;
            app.Spikes.amplitudePP = s.amplitudePP;
            app.Spikes.clusterIds = s.clusterIds;
            app.Classification = s.classification;
            app.refreshClassTable();
            app.updateClusterList();
            app.setSelectedClusters(intersect(s.selection, app.Spikes.clusterIds));
            app.refreshAll();
        end

        function undo(app)
            if isempty(app.UndoStack)
                return
            end
            app.RedoStack{end + 1} = app.snapshotState();
            app.restoreState(app.UndoStack{end});
            app.UndoStack(end) = [];
            app.updateUndoMenus();
            app.setInfo("Undo.");
        end

        function redo(app)
            if isempty(app.RedoStack)
                return
            end
            app.UndoStack{end + 1} = app.snapshotState();
            app.restoreState(app.RedoStack{end});
            app.RedoStack(end) = [];
            app.updateUndoMenus();
            app.setInfo("Redo.");
        end

        function updateUndoMenus(app)
            app.UndoMenu.Enable = matlab.lang.OnOffSwitchState(~isempty(app.UndoStack));
            app.RedoMenu.Enable = matlab.lang.OnOffSwitchState(~isempty(app.RedoStack));
        end

        function reloadData(app)
            if isfield(app.Spikes, "phyDir") && ~isempty(app.Spikes.phyDir)
                app.loadPhyFolder(app.Spikes.phyDir);
            else
                app.setInfo("Nothing to reload.");
            end
        end

        function buildControls(app, outer)
            panel = uigridlayout(outer, [7 1]);
            panel.Layout.Row = 1;
            panel.Layout.Column = 1;
            panel.RowHeight = {"fit", "1x", "fit", "fit", "fit", "1x", "fit"};
            panel.Padding = [6 6 6 6];
            panel.RowSpacing = 8;

            app.buildFileSection(panel, 1);
            app.buildClusterSection(panel, 2);
            app.buildLabelSection(panel, 3);
            app.buildCurateSection(panel, 4);
            app.buildThresholdSection(panel, 5);
            app.buildTableSection(panel, 6);
            app.buildSaveSection(panel, 7);
        end

        function buildFileSection(app, panel, row)
            g = uigridlayout(panel, [1 3]);
            g.Layout.Row = row;
            g.Padding = [0 0 0 0];
            g.ColumnWidth = {"1x", "1x", 70};
            uibutton(g, Text="Load Phy folder...", ...
                ButtonPushedFcn=@(~, ~) app.promptForSource());
            uibutton(g, Text="Load file (.mat)...", ...
                ButtonPushedFcn=@(~, ~) app.loadFileDialog());
            uibutton(g, Text="Reload", ...
                ButtonPushedFcn=@(~, ~) app.reloadData());
        end

        function buildClusterSection(app, panel, row)
            g = uigridlayout(panel, [2 1]);
            g.Layout.Row = row;
            g.RowHeight = {20, "1x"};
            g.Padding = [0 0 0 0];
            g.RowSpacing = 2;
            header = uigridlayout(g, [1 3]);
            header.Padding = [0 0 0 0];
            header.ColumnWidth = {"fit", "fit", "1x"};
            uilabel(header, Text="Clusters", FontWeight="bold");
            uilabel(header, Text="sort:", HorizontalAlignment="right");
            uidropdown(header, Items=["ID", "n", "Label"], ...
                ItemsData=["id", "count", "label"], Value="id", ...
                ValueChangedFcn=@(s, ~) app.setClusterSort(s.Value));
            app.ClusterTable = uitable(g, ...
                ColumnName={'', 'ID', 'n', 'Label'}, ...
                ColumnWidth={24, 34, 55, 'auto'}, ...
                RowName={}, Multiselect="on", SelectionType="row", ...
                SelectionChangedFcn=@(s, e) app.onClusterChanged(s, e));
        end

        function buildLabelSection(app, panel, row)
            g = uigridlayout(panel, [2 4]);
            g.Layout.Row = row;
            g.RowHeight = {16, "fit"};
            g.Padding = [0 0 0 0];
            g.RowSpacing = 2;
            title = uilabel(g, Text="Label selected", FontWeight="bold");
            title.Layout.Row = 1; title.Layout.Column = [1 4];
            labels = ["SU", "MU", "Noise", "Unsorted"];
            for k = 1:numel(labels)
                b = uibutton(g, Text=labels(k), ...
                    ButtonPushedFcn=@(~, ~) app.classifySelected(labels(k)));
                b.Layout.Row = 2; b.Layout.Column = k;
            end
        end

        function buildCurateSection(app, panel, row)
            g = uigridlayout(panel, [2 2]);
            g.Layout.Row = row;
            g.RowHeight = {"fit", "fit"};
            g.Padding = [0 0 0 0];
            g.RowSpacing = 4;
            curate = uibutton(g, Text="Curation", BackgroundColor=[0.85 0.9 1], ...
                ButtonPushedFcn=@(~, ~) app.launchCurateTool());
            curate.Layout.Row = 1; curate.Layout.Column = [1 2];
            m = uibutton(g, Text="Merge selected", ...
                ButtonPushedFcn=@(~, ~) app.mergeSelected());
            m.Layout.Row = 2; m.Layout.Column = 1;
            r = uibutton(g, Text="Realign selected", ...
                ButtonPushedFcn=@(~, ~) app.realignSelected());
            r.Layout.Row = 2; r.Layout.Column = 2;
        end

        function buildThresholdSection(app, panel, row)
            g = uigridlayout(panel, [2 2]);
            g.Layout.Row = row;
            g.RowHeight = {16, "fit"};
            g.Padding = [0 0 0 0];
            g.RowSpacing = 2;
            title = uilabel(g, Text="Trace threshold (waveform crossing)", ...
                FontWeight="bold");
            title.Layout.Row = 1; title.Layout.Column = [1 2];
            b1 = uibutton(g, Text="Add threshold", ...
                ButtonPushedFcn=@(~, ~) app.addThresholdLine());
            b1.Layout.Row = 2; b1.Layout.Column = 1;
            b2 = uibutton(g, Text="Clear thresholds", ...
                ButtonPushedFcn=@(~, ~) app.clearThresholdLines());
            b2.Layout.Row = 2; b2.Layout.Column = 2;
        end

        function buildTableSection(app, panel, row)
            g = uigridlayout(panel, [2 1]);
            g.Layout.Row = row;
            g.RowHeight = {16, "1x"};
            g.Padding = [0 0 0 0];
            g.RowSpacing = 2;
            uilabel(g, Text="Cluster labels", FontWeight="bold");
            app.ClassTable = uitable(g, Data=app.Classification, ...
                ColumnEditable=[false true true], ...
                CellEditCallback=@(~, e) app.onClassEdited(e));
        end

        function buildSaveSection(app, panel, row)
            g = uigridlayout(panel, [1 1]);
            g.Layout.Row = row;
            g.Padding = [0 0 0 0];
            uibutton(g, Text="Save sorting + labels", ...
                ButtonPushedFcn=@(~, ~) app.saveSorting());
        end

        function buildPlots(app, outer)
            grid = uigridlayout(outer, [6 3]);
            grid.Layout.Row = 1;
            grid.Layout.Column = 2;
            % Rows 1-2: top panels; rows 3-4: bottom panels; 5: trace; 6: nav.
            grid.RowHeight = {"1x", "1x", "1x", "1x", "1.4x", 30};

            % Top panels (two rows tall): Mean | ISI | Amplitude.
            app.MeanAxes = uiaxes(grid);
            app.MeanAxes.Layout.Row = [1 2]; app.MeanAxes.Layout.Column = 1;
            app.ISIAxes = uiaxes(grid);
            app.ISIAxes.Layout.Row = [1 2]; app.ISIAxes.Layout.Column = 2;
            app.AmplitudeAxes = uiaxes(grid);
            app.AmplitudeAxes.Layout.Row = [1 2]; app.AmplitudeAxes.Layout.Column = 3;

            % Bottom panels (two rows tall): Waveforms | Features | correlograms.
            app.WaveformAxes = uiaxes(grid);
            app.WaveformAxes.Layout.Row = [3 4]; app.WaveformAxes.Layout.Column = 1;
            app.FeatureAxes = uiaxes(grid);
            app.FeatureAxes.Layout.Row = [3 4]; app.FeatureAxes.Layout.Column = 2;
            app.CorrLayout = uigridlayout(grid, [1 1]);
            app.CorrLayout.Layout.Row = [3 4]; app.CorrLayout.Layout.Column = 3;
            app.CorrLayout.Padding = [2 2 2 2];
            app.CorrLayout.RowSpacing = 2; app.CorrLayout.ColumnSpacing = 2;

            app.TraceAxes = uiaxes(grid);
            app.TraceAxes.Layout.Row = 5; app.TraceAxes.Layout.Column = [1 3];
            addlistener(app.TraceAxes, "XLim", "PostSet", ...
                @(~, ~) app.onTraceXLimChanged());

            % Navigation row: Prev spike | slider | Next spike.
            nav = uigridlayout(grid, [1 3]);
            nav.Layout.Row = 6; nav.Layout.Column = [1 3];
            nav.ColumnWidth = {130, "1x", 130};
            nav.Padding = [0 0 0 0];
            uibutton(nav, Text="< Prev spike", ...
                ButtonPushedFcn=@(~, ~) app.jumpSpike(-1));
            app.TraceSlider = uislider(nav, Limits=[0 1], Value=0.5, ...
                MajorTicks=[], ValueChangedFcn=@(s, e) app.onTraceSlider(s, e));
            uibutton(nav, Text="Next spike >", ...
                ButtonPushedFcn=@(~, ~) app.jumpSpike(1));
        end

        function jumpSpike(app, direction)
            %JUMPSPIKE Centre the trace on the prev/next spike of selected clusters.
            if ~app.ensureTraces() || isnan(app.TraceCenter)
                return
            end
            st = sort(app.Spikes.spikeTimes( ...
                ismember(app.Spikes.clusters, app.selectedClusters())));
            if isempty(st)
                return
            end
            if direction > 0
                target = st(find(st > app.TraceCenter + 1, 1, "first"));
            else
                target = st(find(st < app.TraceCenter - 1, 1, "last"));
            end
            if isempty(target)
                return
            end
            app.TraceCenter = double(target);   % keep the same window length
            app.clampTraceCenter();
            app.updateTrace();
        end

        function loadFileDialog(app)
            [f, d] = uigetfile({'*.mat', "MATLAB spike files"}, ...
                "Select a sorted-spikes .mat file");
            figure(app.UIFigure);
            if isequal(f, 0)
                return
            end
            try
                model = load(fullfile(d, f));
                app.setData(model);
            catch err
                app.setInfo("Could not load as a spike model (needs clusters/" + ...
                    "spikeTimes/waveforms fields): " + err.message);
            end
        end
    end
end

% =========================================================================
% Local helper functions
% =========================================================================
function name = addonDisplayName(file)
    %ADDONDISPLAYNAME Human-readable name for an add-on script.
    %   Uses a "% AddonName: ..." tag if present, else prettifies the file name.
    lines = string(splitlines(fileread(file)));
    tag = regexp(lines, "^\s*%\s*AddonName:\s*(.+)$", "tokens", "once");
    tag = tag(~cellfun(@isempty, tag));
    if ~isempty(tag)
        name = strtrim(string(tag{1}));
        return
    end
    [~, base] = fileparts(file);
    words = split(string(base), "_");
    words(1) = upper(extractBefore(words(1), 2)) + extractAfter(words(1), 1);
    name = join(words, " ");
end

function s = missingToEmpty(s)
    if ismissing(s)
        s = "";
    else
        s = string(s);
    end
end

function t = emptyClassificationTable()
    t = table('Size', [0 3], ...
        'VariableTypes', ["double", "string", "string"], ...
        'VariableNames', ["ClusterID", "Label", "Note"]);
end

function m = normalizeModel(model)
    %NORMALIZEMODEL Coerce a loaded model into the fields the app relies on.
    m = model;
    if ~isfield(m, "sourceType")
        m.sourceType = "mat";
    end
    if ~isfield(m, "name")
        if isfield(m, "phyDir")
            [~, parent] = fileparts(fileparts(m.phyDir));
            [~, leaf] = fileparts(m.phyDir);
            m.name = string(parent) + "/" + string(leaf);
        else
            m.name = "spikes";
        end
    end
    m.clusters = double(m.clusters(:));
    m.spikeTimes = double(m.spikeTimes(:));
    m.clusterIds = unique(m.clusters);

    if isfield(m, "waveforms") && ~isempty(m.waveforms)
        m.wfMinADC = double(min(m.waveforms, [], 2));
        m.wfMaxADC = double(max(m.waveforms, [], 2));
        m.amplitudePP = m.wfMaxADC - m.wfMinADC;
    else
        [m.wfMinADC, m.wfMaxADC] = deal(zeros(numel(m.clusters), 1));
        if isfield(m, "amplitudes")
            m.amplitudePP = double(m.amplitudes(:));
        else
            m.amplitudePP = zeros(numel(m.clusters), 1);
        end
    end

    % Microvolts per ADC count. These exports declare gainToUV as ADC counts
    % per microvolt (dividing yields physiological amplitudes), so invert it.
    if isfield(m, "gainToUV") && isfinite(m.gainToUV) && m.gainToUV > 0
        m.uvPerADC = 1 / m.gainToUV;
    else
        m.uvPerADC = 1;
    end

    if ~isfield(m, "waveformWindow")
        m.waveformWindow = [floor(size(m.waveforms, 2) / 2), ...
            ceil(size(m.waveforms, 2) / 2) - 1];
    end
    if ~isfield(m, "clusterLabels")
        m.clusterLabels = repmat("", numel(m.clusterIds), 1);
        if isfield(m, "clusterTable") && ismember("group", ...
                string(m.clusterTable.Properties.VariableNames))
            for k = 1:numel(m.clusterIds)
                row = m.clusterTable.cluster_id == m.clusterIds(k);
                if any(row)
                    m.clusterLabels(k) = string(m.clusterTable.group(find(row, 1)));
                end
            end
        end
    end
    if ~isfield(m, "numSamples"), m.numSamples = NaN; end
    if ~isfield(m, "numChannels"), m.numChannels = 1; end
    if ~isfield(m, "dataChannel"), m.dataChannel = 0; end
    if ~isfield(m, "samplingRate"), m.samplingRate = 30000; end
end
