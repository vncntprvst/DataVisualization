classdef PETHTool < handle
    %PETHTOOL Peri-event time histograms for the selected clusters.
    %   PETHTool(app) opens a companion window (main app -> Tools) that aligns the
    %   selected clusters' spikes to behavioural events and shows a histogram PETH
    %   per cluster x event.
    %
    %   Events are loaded from events.mat (or events.tsv/.csv) in the Phy folder:
    %   each event name maps to a vector of alignment times in recording seconds
    %   (see LoadEvents). Pick which events to show, the bin size and the
    %   pre/post window; optionally split an event's PETH by its condition label
    %   (e.g. saccade direction).
    %
    %   See also LoadEvents, SpikeVisualizationApp.

    properties (Access = private)
        App
        Events struct = struct()   % from LoadEvents ([] fields if none loaded)

        UIFigure               matlab.ui.Figure
        EventList              matlab.ui.control.ListBox
        BinField               matlab.ui.control.NumericEditField
        PreField               matlab.ui.control.NumericEditField
        PostField              matlab.ui.control.NumericEditField
        SplitBox               matlab.ui.control.CheckBox
        PlotLayout             matlab.ui.container.GridLayout
        InfoLabel              matlab.ui.control.Label
        LastField              % scratch handle for addNumeric
    end

    properties (Constant, Access = private)
        MaxClusters = 6
        MaxEvents = 6
    end

    methods
        function tool = PETHTool(app)
            arguments
                app (1,1) SpikeVisualizationApp
            end
            tool.App = app;
            tool.buildUI();
            tool.tryLoadEvents(app.Spikes.phyDir);
            tool.refresh();
        end

        function delete(tool)
            if ~isempty(tool.UIFigure) && isvalid(tool.UIFigure)
                delete(tool.UIFigure);
            end
        end
    end

    methods (Access = private)
        function buildUI(tool)
            tool.UIFigure = uifigure(Name="PETH / event alignment", ...
                Position=[120 120 1180 700], ...
                CloseRequestFcn=@(~, ~) delete(tool));
            outer = uigridlayout(tool.UIFigure, [2 2]);
            outer.RowHeight = {"1x", 20};
            outer.ColumnWidth = {250, "1x"};

            ctl = uigridlayout(outer, [13 2]);
            ctl.Layout.Row = 1; ctl.Layout.Column = 1;
            ctl.RowHeight = {26, 16, "1x", 16, 26, 16, 26, 16, 26, 26, 26, 26, "fit"};
            ctl.Padding = [6 6 6 6];

            b = uibutton(ctl, Text="Load events", ...
                ButtonPushedFcn=@(~, ~) tool.loadDialog());
            b.Layout.Row = 1; b.Layout.Column = [1 2];
            lbl = uilabel(ctl, Text="Events", FontWeight="bold");
            lbl.Layout.Row = 2; lbl.Layout.Column = [1 2];
            tool.EventList = uilistbox(ctl, Multiselect="on", Items={'(none)'});
            tool.EventList.Layout.Row = 3; tool.EventList.Layout.Column = [1 2];

            tool.addNumeric(ctl, 4, 5, "Bin (ms)", 20);
            tool.BinField = tool.LastField;
            tool.addNumeric(ctl, 6, 7, "Pre (ms)", 500);
            tool.PreField = tool.LastField;
            tool.addNumeric(ctl, 8, 9, "Post (ms)", 500);
            tool.PostField = tool.LastField;

            tool.SplitBox = uicheckbox(ctl, Text="Split by condition");
            tool.SplitBox.Layout.Row = 10; tool.SplitBox.Layout.Column = [1 2];
            rb = uibutton(ctl, Text="Refresh (from app selection)", ...
                BackgroundColor=[0.85 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.refresh());
            rb.Layout.Row = 11; rb.Layout.Column = [1 2];

            tool.PlotLayout = uigridlayout(outer, [1 1]);
            tool.PlotLayout.Layout.Row = 1; tool.PlotLayout.Layout.Column = 2;
            tool.PlotLayout.Padding = [2 2 2 2];

            tool.InfoLabel = uilabel(outer, Text="");
            tool.InfoLabel.Layout.Row = 2; tool.InfoLabel.Layout.Column = [1 2];
        end

        function addNumeric(tool, parent, labelRow, fieldRow, text, value)
            l = uilabel(parent, Text=text);
            l.Layout.Row = labelRow; l.Layout.Column = [1 2];
            f = uieditfield(parent, "numeric", Value=value, ...
                ValueChangedFcn=@(~, ~) tool.refresh());
            f.Layout.Row = fieldRow; f.Layout.Column = [1 2];
            tool.LastField = f;
        end

        function tryLoadEvents(tool, source)
            try
                tool.setEvents(LoadEvents(string(source)));
            catch err
                tool.setInfo("No events loaded: " + err.message + ...
                    "  Use 'Load events'.");
            end
        end

        function loadDialog(tool)
            [f, d] = uigetfile({'*.mat;*.tsv;*.csv', "Event files"}, ...
                "Select an events file");
            figure(tool.UIFigure);
            if isequal(f, 0)
                return
            end
            try
                tool.setEvents(LoadEvents(string(fullfile(d, f))));
                tool.refresh();
            catch err
                tool.setInfo("Load failed: " + err.message);
            end
        end

        function setEvents(tool, ev)
            tool.Events = ev;
            names = ev.names;
            tool.EventList.Items = cellstr(names);
            tool.EventList.ItemsData = names;
            wanted = intersect(["target", "saccade", "reward"], names, "stable");
            if isempty(wanted)
                wanted = names(1:min(3, numel(names)));
            end
            tool.EventList.Value = wanted;
            tool.setInfo("Loaded " + numel(names) + " events from " + ev.source);
        end

        function refresh(tool)
            delete(tool.PlotLayout.Children);
            if ~isfield(tool.Events, "events")
                tool.setInfo("Load an events file first.");
                return
            end
            clusters = tool.App.selectedClusterIds();
            clusters = clusters(1:min(numel(clusters), tool.MaxClusters));
            events = string(tool.EventList.Value);
            events = events(1:min(numel(events), tool.MaxEvents));
            if isempty(clusters) || isempty(events)
                tool.setInfo("Select cluster(s) in the app and event(s) here.");
                return
            end

            nr = numel(clusters);
            nc = numel(events);
            tool.PlotLayout.RowHeight = repmat({"1x"}, 1, nr);
            tool.PlotLayout.ColumnWidth = repmat({"1x"}, 1, nc);
            fs = tool.App.Spikes.samplingRate;
            preS = tool.PreField.Value / 1000;
            postS = tool.PostField.Value / 1000;
            binS = max(tool.BinField.Value, 1) / 1000;

            for r = 1:nr
                spikeSec = tool.App.Spikes.spikeTimes( ...
                    tool.App.Spikes.clusters == clusters(r)) / fs;
                for c = 1:nc
                    ax = uiaxes(tool.PlotLayout);
                    ax.Layout.Row = r; ax.Layout.Column = c;
                    tool.plotCell(ax, spikeSec, clusters(r), events(c), ...
                        preS, postS, binS, r, c, nr);
                end
            end
            tool.setInfo(sprintf("PETH: %d cluster(s) x %d event(s), " + ...
                "%g ms bins, [-%g +%g] ms.", nr, nc, tool.BinField.Value, ...
                tool.PreField.Value, tool.PostField.Value));
        end

        function plotCell(tool, ax, spikeSec, cid, evName, preS, postS, binS, r, c, nr)
            fieldName = matlab.lang.makeValidName(evName);
            align = tool.Events.events.(fieldName);
            edges = -preS:binS:postS;
            centres = (edges(1:end-1) + binS / 2) * 1000;   % ms

            splitOn = tool.SplitBox.Value && isfield(tool.Events.cond, fieldName) ...
                && ~isempty(tool.Events.cond.(fieldName));
            hold(ax, "on");
            if splitOn
                cond = tool.Events.cond.(fieldName);
                groups = unique(cond);
                pal = lines(max(numel(groups), 3));
                for gi = 1:numel(groups)
                    [meanR, sem] = condSDF(spikeSec, align(cond == groups(gi)), ...
                        edges, binS);
                    band = patch(ax, [centres fliplr(centres)], ...
                        [meanR - sem fliplr(meanR + sem)], pal(gi, :), ...
                        EdgeColor="none", FaceAlpha=0.15);
                    band.Annotation.LegendInformation.IconDisplayStyle = "off";
                    plot(ax, centres, meanR, Color=pal(gi, :), LineWidth=1.5, ...
                        DisplayName=groups(gi));
                end
                if r == 1
                    legend(ax, Location="best", Box="off");
                end
            else
                rate = pethRate(spikeSec, align, edges, binS);
                bar(ax, centres, rate, 1, FaceColor=tool.App.clusterColor(cid), ...
                    EdgeColor="none");
            end
            xl = xline(ax, 0, "-", Color=[0.4 0.4 0.4]);
            xl.Annotation.LegendInformation.IconDisplayStyle = "off";   % keep out of legend
            hold(ax, "off");
            xlim(ax, [-preS postS] * 1000);
            if r == 1
                title(ax, evName + "  (n=" + numel(align) + ")", Interpreter="none");
            end
            if c == 1
                ylabel(ax, sprintf("c%d  (Hz)", cid));
            end
            if r == nr
                xlabel(ax, "Time (ms)");
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
function rate = pethRate(spikeSec, alignSec, edges, binS)
    %PETHRATE Mean firing rate (Hz) in each bin, averaged over alignment times.
    counts = zeros(1, numel(edges) - 1);
    span = [edges(1), edges(end)];
    for i = 1:numel(alignSec)
        d = spikeSec - alignSec(i);
        d = d(d >= span(1) & d <= span(2));
        counts = counts + histcounts(d, edges);
    end
    rate = counts / max(numel(alignSec), 1) / binS;
end

function [meanR, sem] = condSDF(spikeSec, alignSec, edges, binS)
    %CONDSDF Smoothed density (mean +/- SEM across trials) in Hz.
    nBins = numel(edges) - 1;
    if isempty(alignSec)
        [meanR, sem] = deal(zeros(1, nBins));
        return
    end
    perTrial = zeros(numel(alignSec), nBins);
    span = [edges(1), edges(end)];
    for i = 1:numel(alignSec)
        d = spikeSec - alignSec(i);
        perTrial(i, :) = histcounts(d(d >= span(1) & d <= span(2)), edges);
    end
    sdf = smoothRows(perTrial / binS);   % Hz, Gaussian-smoothed per trial
    meanR = mean(sdf, 1);
    if size(sdf, 1) > 1
        sem = std(sdf, 0, 1) / sqrt(size(sdf, 1));
    else
        sem = zeros(1, nBins);
    end
end

function y = smoothRows(m)
    %SMOOTHROWS Gaussian smoothing (sigma 2 bins) with replicate padding, per row.
    g = exp(-0.5 * ((-6:6) / 2) .^ 2);
    g = g / sum(g);
    pad = 6;
    y = zeros(size(m));
    for r = 1:size(m, 1)
        x = m(r, :);
        xp = [repmat(x(1), 1, pad), x, repmat(x(end), 1, pad)];
        yc = conv(xp, g, "same");
        y(r, :) = yc(pad + 1:end - pad);
    end
end
