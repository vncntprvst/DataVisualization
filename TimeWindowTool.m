classdef TimeWindowTool < handle
    %TIMEWINDOWTOOL Restrict curation to a draggable time window.
    %   TimeWindowTool(app, curate) opens a companion window (launched from the
    %   Curation widget's Tools menu) for splitting units that swap amplitude
    %   rank over time on a single electrode.
    %
    %   The top panel is amplitude (peak-to-peak) vs time for the curated
    %   clusters, with two draggable vertical bars marking the window start and
    %   stop. The bottom panel is a raw-trace excerpt at the bar you last moved,
    %   so you can verify the cut. Moving a bar restricts the Curation widget to
    %   the clusters' spikes within [start, stop] — isolate a band where the two
    %   amplitudes are separated, split there, then step to the next band.
    %
    %   Buttons snap the start to 0, the stop to the end, or revert to the whole
    %   recording.
    %
    %   See also CurationTool, SpikeVisualizationApp.

    properties (Access = private)
        App
        Curate                 % CurationTool
        ClusterIds double

        UIFigure               matlab.ui.Figure
        AmpAxes                matlab.ui.control.UIAxes
        TraceAxes              matlab.ui.control.UIAxes
        InfoLabel              matlab.ui.control.Label
        StartLine              images.roi.Line
        StopLine               images.roi.Line

        Suppress logical = false
        FocusSec double        % time the trace excerpt is centred on
        AmpYLim double

        CurationSelIdx double = []   % app spike indices selected in Curation
        PlotIdx double = []          % app spike index behind each amp point
        PlotX double = []            % time (s) of each amp point
        PlotY double = []            % peak-to-peak (uV) of each amp point
    end

    methods
        function tool = TimeWindowTool(app, curate)
            arguments
                app (1,1) SpikeVisualizationApp
                curate (1,1) CurationTool
            end
            tool.App = app;
            tool.Curate = curate;
            tool.ClusterIds = curate.shownClusters();
            tool.CurationSelIdx = curate.globalSelection();
            w = app.timeWindowSec();
            tool.FocusSec = w(1);
            tool.buildUI();
            tool.redrawAmp();
            tool.plotTrace();
            tool.updateInfo();
        end

        function delete(tool)
            if ~isempty(tool.UIFigure) && isvalid(tool.UIFigure)
                delete(tool.UIFigure);
            end
        end

        function h = figureHandle(tool)
            h = tool.UIFigure;
        end

        function setCurationSelection(tool, globalIdx)
            %SETCURATIONSELECTION Reflect the Curation selection here (live link).
            tool.CurationSelIdx = globalIdx;
            tool.redrawAmp();
        end
    end

    methods (Access = private)
        function dur = duration(tool)
            dur = tool.App.Spikes.numSamples / tool.App.Spikes.samplingRate;
        end

        function buildUI(tool)
            tool.UIFigure = uifigure(Name="Time-window split", ...
                Position=[160 120 1000 620], ...
                CloseRequestFcn=@(~, ~) delete(tool));
            g = uigridlayout(tool.UIFigure, [3 1]);
            g.RowHeight = {"1x", "1x", 34};
            g.Padding = [6 6 6 6];

            tool.AmpAxes = uiaxes(g);
            tool.AmpAxes.Layout.Row = 1;
            tool.TraceAxes = uiaxes(g);
            tool.TraceAxes.Layout.Row = 2;

            btns = uigridlayout(g, [1 5]);
            btns.Layout.Row = 3;
            btns.Padding = [0 0 0 0];
            btns.ColumnWidth = {130, 120, 120, 160, "1x"};
            uibutton(btns, Text="Lasso select", BackgroundColor=[0.9 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.lassoSelect());
            uibutton(btns, Text="Snap start -> 0", ...
                ButtonPushedFcn=@(~, ~) tool.snapStart());
            uibutton(btns, Text="Snap stop -> end", ...
                ButtonPushedFcn=@(~, ~) tool.snapStop());
            uibutton(btns, Text="Revert (whole recording)", ...
                ButtonPushedFcn=@(~, ~) tool.revert());
            tool.InfoLabel = uilabel(btns, Text="");
        end

        function redrawAmp(tool)
            % Redraw the amplitude plot then the boundary bars (cla drops them).
            tool.plotAmplitude();
            tool.drawBoundaries();
        end

        function plotAmplitude(tool)
            ax = tool.AmpAxes;
            cla(ax);
            hold(ax, "on");
            s = tool.App.Spikes;
            uv = 1;
            if isfield(s, "uvPerADC")
                uv = s.uvPerADC;
            end
            allIdx = find(ismember(s.clusters, tool.ClusterIds));
            tool.PlotIdx = allIdx;
            tool.PlotX = s.spikeTimes(allIdx) / s.samplingRate;
            tool.PlotY = s.amplitudePP(allIdx) * uv;
            for cid = tool.ClusterIds
                m = s.clusters(allIdx) == cid;
                scatter(ax, tool.PlotX(m), tool.PlotY(m), 4, ...
                    tool.App.clusterColor(cid), "filled", MarkerFaceAlpha=0.3);
            end
            hi = ismember(allIdx, tool.CurationSelIdx);   % Curation selection
            if any(hi)
                scatter(ax, tool.PlotX(hi), tool.PlotY(hi), 16, [0 0 0], ...
                    "o", LineWidth=0.75);
            end
            hold(ax, "off");
            title(ax, sprintf("Amplitude vs time (drag bars; %d selected)", sum(hi)));
            xlabel(ax, "Time (s)");
            ylabel(ax, "Peak-to-peak (\muV)");
            xlim(ax, [0 tool.duration()]);
            tool.AmpYLim = ylim(ax);
        end

        function lassoSelect(tool)
            %LASSOSELECT Freehand-select amplitude points -> Curation selection.
            roi = drawfreehand(tool.AmpAxes, Color=[0.1 0.1 0.1]);
            if isempty(roi) || ~isvalid(roi) || size(roi.Position, 1) < 3
                return
            end
            poly = roi.Position;
            delete(roi);
            inside = inpolygon(tool.PlotX, tool.PlotY, poly(:, 1), poly(:, 2));
            if isvalid(tool.Curate)
                tool.Curate.selectGlobal(tool.PlotIdx(inside));
            end
        end

        function drawBoundaries(tool)
            tool.Suppress = true;
            delete(tool.StartLine(isvalid(tool.StartLine)));
            delete(tool.StopLine(isvalid(tool.StopLine)));
            w = tool.App.timeWindowSec();
            yl = tool.AmpYLim;
            tool.StartLine = images.roi.Line(tool.AmpAxes, ...
                Position=[w(1) yl(1); w(1) yl(2)], Color=[0 0.6 0], LineWidth=1.5);
            tool.StopLine = images.roi.Line(tool.AmpAxes, ...
                Position=[w(2) yl(1); w(2) yl(2)], Color=[0.85 0 0], LineWidth=1.5);
            addlistener(tool.StartLine, "MovingROI", ...
                @(s, ~) tool.constrainVertical(s));
            addlistener(tool.StopLine, "MovingROI", ...
                @(s, ~) tool.constrainVertical(s));
            addlistener(tool.StartLine, "ROIMoved", @(~, ~) tool.onMoved("start"));
            addlistener(tool.StopLine, "ROIMoved", @(~, ~) tool.onMoved("stop"));
            tool.Suppress = false;
        end

        function constrainVertical(tool, roi)
            if tool.Suppress
                return
            end
            x = min(max(mean(roi.Position(:, 1)), 0), tool.duration());
            roi.Position = [x tool.AmpYLim(1); x tool.AmpYLim(2)];
        end

        function onMoved(tool, which)
            if tool.Suppress
                return
            end
            tool.FocusSec = mean(tool.movedLine(which).Position(:, 1));
            tool.applyWindow();
            tool.plotTrace();
        end

        function roi = movedLine(tool, which)
            if which == "start"
                roi = tool.StartLine;
            else
                roi = tool.StopLine;
            end
        end

        function applyWindow(tool)
            a = mean(tool.StartLine.Position(:, 1));
            b = mean(tool.StopLine.Position(:, 1));
            tool.App.setTimeWindow(sort([a b]));   % refreshes the Curation tool
            tool.updateInfo();
        end

        function plotTrace(tool)
            ax = tool.TraceAxes;
            cla(ax);
            [tSec, uv] = tool.App.traceExcerpt(tool.FocusSec, 0.5);
            if isempty(tSec)
                title(ax, "Raw trace unavailable");
                return
            end
            plot(ax, tSec, uv, Color=[0 0 0 0.8], LineWidth=0.1);
            hold(ax, "on");
            xline(ax, tool.FocusSec, "-", Color=[0.85 0 0], LineWidth=1.5);
            s = tool.App.Spikes;
            yb = double(min(uv));
            for cid = tool.ClusterIds
                st = s.spikeTimes(s.clusters == cid) / s.samplingRate;
                st = st(st >= tSec(1) & st <= tSec(end));
                plot(ax, st, repmat(yb, size(st)), "^", ...
                    Color=tool.App.clusterColor(cid), ...
                    MarkerFaceColor=tool.App.clusterColor(cid), ...
                    MarkerSize=4, LineStyle="none");
            end
            hold(ax, "off");
            title(ax, sprintf("Raw trace at boundary (%.2f s)", tool.FocusSec));
            xlabel(ax, "Time (s)");
            ylabel(ax, "Voltage (\muV)");
            axis(ax, "tight");
        end

        function snapStart(tool)
            tool.Suppress = true;
            tool.StartLine.Position = [0 tool.AmpYLim(1); 0 tool.AmpYLim(2)];
            tool.Suppress = false;
            tool.FocusSec = 0;
            tool.applyWindow();
            tool.plotTrace();
        end

        function snapStop(tool)
            d = tool.duration();
            tool.Suppress = true;
            tool.StopLine.Position = [d tool.AmpYLim(1); d tool.AmpYLim(2)];
            tool.Suppress = false;
            tool.FocusSec = d;
            tool.applyWindow();
            tool.plotTrace();
        end

        function revert(tool)
            tool.App.setTimeWindow([]);   % whole recording
            tool.drawBoundaries();
            tool.updateInfo();
        end

        function updateInfo(tool)
            w = tool.App.timeWindowSec();
            s = tool.App.Spikes;
            tSec = s.spikeTimes / s.samplingRate;
            inWin = ismember(s.clusters, tool.ClusterIds) & ...
                tSec >= w(1) & tSec <= w(2);
            tool.InfoLabel.Text = sprintf("Window %.1f - %.1f s  |  %d spikes", ...
                w(1), w(2), sum(inWin));
        end
    end
end
