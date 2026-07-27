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
    %   The window shows, for the selected cluster(s): overlaid and mean
    %   waveforms, the inter-spike-interval histogram, the autocorrelogram,
    %   amplitude drift over time, a PCA feature scatter, and a scrollable raw
    %   data excerpt with spike markers. Clusters can be labelled (single unit,
    %   multi unit, noise) and the labels saved. The "ISI split..." button opens
    %   ISISplitTool to lasso refractory-violating spikes and split them out.
    %
    %   This is the modern uifigure replacement for the legacy GUIDE
    %   SpikeVisualizationGUI.
    %
    %   See also LoadSpikesPhy, ISISplitTool.

    properties (SetAccess = private)
        Spikes struct = struct()   % normalized data model (see LoadSpikesPhy)
        Classification table       % one row per labelled cluster
    end

    properties (Access = private)
        UIFigure               matlab.ui.Figure
        ClusterList            matlab.ui.control.ListBox
        ClassTable             matlab.ui.control.Table
        WaveformAxes           matlab.ui.control.UIAxes
        MeanAxes               matlab.ui.control.UIAxes
        AmplitudeAxes          matlab.ui.control.UIAxes
        ISIAxes                matlab.ui.control.UIAxes
        ACGAxes                matlab.ui.control.UIAxes
        FeatureAxes            matlab.ui.control.UIAxes
        TraceAxes              matlab.ui.control.UIAxes
        TraceSlider            matlab.ui.control.Slider
        InfoLabel              matlab.ui.control.Label
        SplitTool                          % ISISplitTool handle (companion)

        Traces = []            % memmapfile of the raw data, lazily opened
        TraceCenter double = NaN
        TraceHalfWidth double = NaN
    end

    properties (Constant, Access = private)
        MaxWaveformsToPlot = 200
        ISIMaxMs = 50          % ISI histogram / ACG span in milliseconds
        ISIBinMs = 0.5
        RefractoryMs = 2       % shaded refractory region in the ISI plot
    end

    methods
        function app = SpikeVisualizationApp(source)
            arguments
                source = []    % [] -> ask; folder string -> Phy; struct -> model
            end
            app.Classification = emptyClassificationTable();
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
            if ~isempty(app.SplitTool) && isvalid(app.SplitTool)
                delete(app.SplitTool);
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
        end

        function setData(app, model)
            %SETDATA Adopt a normalized spike model and refresh every view.
            arguments
                app
                model (1,1) struct
            end
            app.Spikes = normalizeModel(model);
            app.Traces = [];   % force re-open of the raw data on next trace draw
            app.TraceCenter = NaN;
            app.Classification = emptyClassificationTable();
            app.seedClassificationFromModel();
            app.updateClusterList();
            % Select the most numerous cluster by default.
            ids = app.Spikes.clusterIds;
            counts = arrayfun(@(c) sum(app.Spikes.clusters == c), ids);
            [~, best] = max(counts);
            app.ClusterList.Value = ids(best);
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
            app.plotWaveforms();
            app.plotMean();
            app.plotAmplitude();
            app.plotISI();
            app.plotACG();
            app.plotFeatures();
            app.plotTrace();
        end

        function cids = selectedClusters(app)
            cids = app.ClusterList.Value;
            if isempty(cids)
                cids = double.empty(0, 1);
            else
                cids = double(cids(:));
            end
        end

        function idx = spikesOf(app, cids)
            idx = ismember(app.Spikes.clusters, cids);
        end

        function plotWaveforms(app)
            ax = app.WaveformAxes;
            cla(ax);
            hold(ax, "on");
            wf = app.Spikes.waveforms;
            t = app.waveformTimeMs();
            for cid = app.selectedClusters()'
                idx = find(app.Spikes.clusters == cid);
                if isempty(idx) || isempty(wf)
                    continue
                end
                idx = idx(round(linspace(1, numel(idx), ...
                    min(numel(idx), app.MaxWaveformsToPlot))));
                plot(ax, t, wf(idx, :)', Color=[app.clusterColor(cid) 0.15]);
            end
            hold(ax, "off");
            title(ax, "Waveforms");
            xlabel(ax, "Time (ms)");
            ylabel(ax, "ADC");
            axis(ax, "tight");
        end

        function plotMean(app)
            ax = app.MeanAxes;
            cla(ax);
            hold(ax, "on");
            wf = app.Spikes.waveforms;
            t = app.waveformTimeMs();
            for cid = app.selectedClusters()'
                idx = app.Spikes.clusters == cid;
                if ~any(idx) || isempty(wf)
                    continue
                end
                m = mean(wf(idx, :), 1);
                sem = std(wf(idx, :), 0, 1) ./ sqrt(sum(idx)) * 1.96;
                c = app.clusterColor(cid);
                fill(ax, [t fliplr(t)], [m - sem fliplr(m + sem)], c, ...
                    EdgeColor="none", FaceAlpha=0.2);
                plot(ax, t, m, Color=c, LineWidth=1.5);
            end
            hold(ax, "off");
            title(ax, "Mean waveform (\pm95% CI)");
            xlabel(ax, "Time (ms)");
            ylabel(ax, "ADC");
            axis(ax, "tight");
        end

        function plotAmplitude(app)
            ax = app.AmplitudeAxes;
            cla(ax);
            hold(ax, "on");
            fs = app.Spikes.samplingRate;
            for cid = app.selectedClusters()'
                idx = app.Spikes.clusters == cid;
                if ~any(idx)
                    continue
                end
                scatter(ax, app.Spikes.spikeTimes(idx) / fs, ...
                    app.Spikes.amplitudePP(idx), 4, app.clusterColor(cid), ...
                    "filled", MarkerFaceAlpha=0.3);
            end
            hold(ax, "off");
            title(ax, "Amplitude vs time (drift)");
            xlabel(ax, "Time (s)");
            ylabel(ax, "Peak-to-peak");
            axis(ax, "tight");
        end

        function plotISI(app)
            ax = app.ISIAxes;
            cla(ax);
            cids = app.selectedClusters();
            if isempty(cids)
                return
            end
            % ISI is only meaningful per cluster; use the most numerous one.
            cid = app.dominantCluster(cids);
            isiMs = app.isiMsOf(cid);
            edges = 0:app.ISIBinMs:app.ISIMaxMs;
            histogram(ax, isiMs, edges, FaceColor=app.clusterColor(cid), ...
                EdgeColor="none");
            xline(ax, app.RefractoryMs, "--", Color=[0.6 0 0]);
            nViol = sum(isiMs < app.RefractoryMs);
            pctViol = 100 * nViol / max(1, numel(isiMs));
            title(ax, sprintf("ISI  (cluster %d)  %.2f%% < %g ms", ...
                cid, pctViol, app.RefractoryMs));
            xlabel(ax, "ISI (ms)");
            ylabel(ax, "Count");
            xlim(ax, [0 app.ISIMaxMs]);
        end

        function plotACG(app)
            ax = app.ACGAxes;
            cla(ax);
            cids = app.selectedClusters();
            if isempty(cids)
                return
            end
            cid = app.dominantCluster(cids);
            [counts, lags] = app.correlogram(cid, cid);
            bar(ax, lags, counts, 1, FaceColor=app.clusterColor(cid), ...
                EdgeColor="none");
            title(ax, sprintf("Autocorrelogram (cluster %d)", cid));
            xlabel(ax, "Lag (ms)");
            ylabel(ax, "Count");
            xlim(ax, [-app.ISIMaxMs app.ISIMaxMs]);
        end

        function plotFeatures(app)
            ax = app.FeatureAxes;
            cla(ax);
            hold(ax, "on");
            for cid = app.selectedClusters()'
                sc = app.featureScores(cid);
                if isempty(sc)
                    continue
                end
                scatter(ax, sc(:, 1), sc(:, 2), 5, app.clusterColor(cid), ...
                    "filled", MarkerFaceAlpha=0.3);
            end
            hold(ax, "off");
            title(ax, "PC features (per-cluster PCA of waveforms)");
            xlabel(ax, "PC1");
            ylabel(ax, "PC2");
            axis(ax, "tight");
        end

        function plotTrace(app)
            ax = app.TraceAxes;
            cla(ax);
            if ~app.ensureTraces()
                title(ax, "Raw trace unavailable (no data file)");
                return
            end
            fs = app.Spikes.samplingRate;
            if isnan(app.TraceCenter)
                app.TraceHalfWidth = round(fs / 2);           % 1 s window
                app.TraceCenter = min(app.Spikes.numSamples / 2, ...
                    app.Spikes.spikeTimes(app.dominantFirstSpike()));
            end
            half = app.TraceHalfWidth;
            first = max(1, round(app.TraceCenter) - half);
            last = min(app.Spikes.numSamples, round(app.TraceCenter) + half);
            ch = app.Spikes.dataChannel + 1;
            excerpt = single(app.Traces.Data.raw(ch, first:last));
            tSec = (first:last) / fs;
            plot(ax, tSec, excerpt, Color=[0 0 0 0.8], LineWidth=0.1);
            hold(ax, "on");
            y = double(min(excerpt));
            for cid = app.selectedClusters()'
                st = app.Spikes.spikeTimes(app.Spikes.clusters == cid);
                st = st(st >= first & st <= last);
                plot(ax, st / fs, repmat(y, size(st)), "^", ...
                    Color=app.clusterColor(cid), MarkerFaceColor=app.clusterColor(cid), ...
                    MarkerSize=4, LineStyle="none");
            end
            hold(ax, "off");
            title(ax, "Raw data excerpt");
            xlabel(ax, "Time (s)");
            axis(ax, "tight");
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
            edges = -app.ISIMaxMs:app.ISIBinMs:app.ISIMaxMs;
            lags = edges(1:end-1) + app.ISIBinMs / 2;
            counts = zeros(1, numel(lags));
            win = app.ISIMaxMs;
            for k = 1:numel(a)
                d = b - a(k);
                d = d(abs(d) <= win);
                counts = counts + histcounts(d, edges);
            end
            if cidA == cidB
                counts(abs(lags) < app.ISIBinMs) = 0;   % remove the zero-lag self peak
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
            ok = true;
        end

        function c = clusterColor(app, cid)
            k = find(app.Spikes.clusterIds == cid, 1);
            if isempty(k)
                k = 1;
            end
            cmap = lines(max(7, numel(app.Spikes.clusterIds)));
            c = cmap(mod(k - 1, size(cmap, 1)) + 1, :);
        end
    end

    % ---------------------------------------------------------- classification
    methods (Access = private)
        function classifySelected(app, label)
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
            app.ClassTable.Data = app.Classification;
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
            %SAVESORTING Write the (possibly edited) cluster assignment back to Phy.
            if ~isfield(app.Spikes, "phyDir") || app.Spikes.sourceType ~= "phy"
                app.setInfo("Save sorting only supported for Phy sources.");
                return
            end
            target = fullfile(app.Spikes.phyDir, "spike_clusters.npy");
            backup = fullfile(app.Spikes.phyDir, "spike_clusters.bak.npy");
            if isfile(target) && ~isfile(backup)
                copyfile(target, backup);
            end
            writeNPY(int32(app.Spikes.clusters), target);
            app.saveClassification();
            app.setInfo("Wrote spike_clusters.npy (backup: spike_clusters.bak.npy).");
        end
    end

    % ---------------------------------------------- public API for split tool
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

        function applySplit(app, spikeIndices, newCluster)
            %APPLYSPLIT Reassign the given global spike indices to a new cluster id.
            %   Called by ISISplitTool once the user confirms a split.
            arguments
                app
                spikeIndices (:, 1) double
                newCluster (1, 1) double
            end
            app.Spikes.clusters(spikeIndices) = newCluster;
            app.Spikes.clusterIds = unique(app.Spikes.clusters);
            app.updateClusterList();
            app.ClusterList.Value = intersect(app.ClusterList.Value, ...
                app.Spikes.clusterIds);
            app.refreshAll();
            app.setInfo(sprintf("Split %d spikes into cluster %d.", ...
                numel(spikeIndices), newCluster));
        end
    end

    % ------------------------------------------------------------- split tool
    methods (Access = private)
        function launchSplitTool(app)
            cids = app.selectedClusters();
            if numel(cids) ~= 1
                app.setInfo("Select exactly one cluster to split.");
                return
            end
            if ~isempty(app.SplitTool) && isvalid(app.SplitTool)
                delete(app.SplitTool);
            end
            app.SplitTool = ISISplitTool(app, cids);
        end
    end

    % --------------------------------------------------------------- plumbing
    methods (Access = private)
        function promptForSource(app)
            folder = uigetdir(pwd, "Select a Phy export folder (contains params.py)");
            figure(app.UIFigure);
            if isequal(folder, 0)
                app.setInfo("No folder selected. Use File controls to load a sorting.");
                return
            end
            app.loadPhyFolder(string(folder));
        end

        function updateClusterList(app)
            ids = app.Spikes.clusterIds;
            items = strings(numel(ids), 1);
            for k = 1:numel(ids)
                n = sum(app.Spikes.clusters == ids(k));
                label = "";
                row = app.Classification.ClusterID == ids(k);
                if any(row)
                    label = "  [" + app.Classification.Label(find(row, 1)) + "]";
                end
                items(k) = sprintf("%d  (n=%d)%s", ids(k), n, label);
            end
            app.ClusterList.Items = cellstr(items);
            app.ClusterList.ItemsData = ids;
        end

        function setInfo(app, msg)
            if ~isempty(app.InfoLabel) && isvalid(app.InfoLabel)
                app.InfoLabel.Text = string(msg);
            end
        end

        function onClusterChanged(app, ~, ~)
            app.TraceCenter = NaN;   % recentre trace on the newly selected cluster
            app.refreshAll();
        end

        function onTraceSlider(app, ~, event)
            if ~app.ensureTraces()
                return
            end
            app.TraceCenter = event.Value * app.Spikes.numSamples;
            app.plotTrace();
        end
    end

    % ----------------------------------------------------------------- buildUI
    methods (Access = private)
        function buildUI(app)
            app.UIFigure = uifigure(Name="SpikeVisualizationApp", ...
                Position=[80 80 1400 860]);
            outer = uigridlayout(app.UIFigure, [2 2]);
            outer.RowHeight = {"1x", 24};
            outer.ColumnWidth = {280, "1x"};

            app.buildControls(outer);
            app.buildPlots(outer);

            app.InfoLabel = uilabel(outer, Text="Ready.", ...
                FontColor=[0.2 0.2 0.2]);
            app.InfoLabel.Layout.Row = 2;
            app.InfoLabel.Layout.Column = [1 2];
        end

        function buildControls(app, outer)
            panel = uigridlayout(outer, [10 2]);
            panel.Layout.Row = 1;
            panel.Layout.Column = 1;
            panel.RowHeight = {28, 28, "1x", 28, 28, 28, 28, 28, "fit", 28};
            panel.ColumnWidth = {"1x", "1x"};

            b1 = uibutton(panel, Text="Load Phy folder...", ...
                ButtonPushedFcn=@(~, ~) app.promptForSource());
            b1.Layout.Row = 1; b1.Layout.Column = [1 2];

            b2 = uibutton(panel, Text="Load file (.mat)...", ...
                ButtonPushedFcn=@(~, ~) app.loadFileDialog());
            b2.Layout.Row = 2; b2.Layout.Column = [1 2];

            app.ClusterList = uilistbox(panel, Multiselect="on", ...
                Items={'(no data)'}, ...
                ValueChangedFcn=@(s, e) app.onClusterChanged(s, e));
            app.ClusterList.Layout.Row = 3;
            app.ClusterList.Layout.Column = [1 2];

            labels = ["Single unit", "SU"; "Multi unit", "MU"; ...
                "Noise", "Noise"; "Unsorted", "Unsorted"];
            for k = 1:size(labels, 1)
                btn = uibutton(panel, Text=labels(k, 1), ...
                    ButtonPushedFcn=@(~, ~) app.classifySelected(labels(k, 2)));
                btn.Layout.Row = 3 + k;
                btn.Layout.Column = mod(k - 1, 2) + 1;
            end

            splitBtn = uibutton(panel, Text="ISI split...", ...
                BackgroundColor=[0.85 0.9 1], ...
                ButtonPushedFcn=@(~, ~) app.launchSplitTool());
            splitBtn.Layout.Row = 6;
            splitBtn.Layout.Column = [1 2];

            app.ClassTable = uitable(panel, Data=app.Classification);
            app.ClassTable.Layout.Row = 9;
            app.ClassTable.Layout.Column = [1 2];

            saveBtn = uibutton(panel, Text="Save sorting + labels", ...
                ButtonPushedFcn=@(~, ~) app.saveSorting());
            saveBtn.Layout.Row = 10;
            saveBtn.Layout.Column = [1 2];
        end

        function buildPlots(app, outer)
            grid = uigridlayout(outer, [4 3]);
            grid.Layout.Row = 1;
            grid.Layout.Column = 2;
            grid.RowHeight = {"1x", "1x", "1x", 28};

            app.WaveformAxes = uiaxes(grid);
            app.WaveformAxes.Layout.Row = 1; app.WaveformAxes.Layout.Column = 1;
            app.MeanAxes = uiaxes(grid);
            app.MeanAxes.Layout.Row = 1; app.MeanAxes.Layout.Column = 2;
            app.AmplitudeAxes = uiaxes(grid);
            app.AmplitudeAxes.Layout.Row = 1; app.AmplitudeAxes.Layout.Column = 3;

            app.ISIAxes = uiaxes(grid);
            app.ISIAxes.Layout.Row = 2; app.ISIAxes.Layout.Column = 1;
            app.ACGAxes = uiaxes(grid);
            app.ACGAxes.Layout.Row = 2; app.ACGAxes.Layout.Column = 2;
            app.FeatureAxes = uiaxes(grid);
            app.FeatureAxes.Layout.Row = 2; app.FeatureAxes.Layout.Column = 3;

            app.TraceAxes = uiaxes(grid);
            app.TraceAxes.Layout.Row = 3; app.TraceAxes.Layout.Column = [1 3];

            app.TraceSlider = uislider(grid, Limits=[0 1], Value=0.5, ...
                MajorTicks=[], ValueChangedFcn=@(s, e) app.onTraceSlider(s, e));
            app.TraceSlider.Layout.Row = 4;
            app.TraceSlider.Layout.Column = [1 3];
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
    if ~isfield(m, "amplitudePP")
        if isfield(m, "waveforms") && ~isempty(m.waveforms)
            m.amplitudePP = double(max(m.waveforms, [], 2) - min(m.waveforms, [], 2));
        elseif isfield(m, "amplitudes")
            m.amplitudePP = double(m.amplitudes(:));
        else
            m.amplitudePP = zeros(numel(m.clusters), 1);
        end
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
