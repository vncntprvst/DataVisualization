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
    %   "Split bands" detects horizontal amplitude bands (drift-tolerant)
    %   among the spikes inside the current [start, stop] window - place the
    %   bars around an epoch to band-split just that stretch. Amplitudes are
    %   divided by the main band's running median over time, and the bands
    %   are the modes of the resulting ratio density. Band 1 is the main
    %   (largest) band.
    %   Pick a band to see its spikes selected in the Curation widget, split
    %   one off with "Band -> new cluster", or give every out-of-band its
    %   own cluster in one go.
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

        BandIdx double = []          % app spike index behind each band label
        BandLabels double = []       % band per BandIdx entry (1 = main band)
        BandDropdown   matlab.ui.control.DropDown
        BandMoveBtn    matlab.ui.control.Button
        BandMoveAllBtn matlab.ui.control.Button
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
                Position=[160 120 1000 660], ...
                CloseRequestFcn=@(~, ~) delete(tool));
            g = uigridlayout(tool.UIFigure, [4 1]);
            g.RowHeight = {"1x", "1x", 34, 34};
            g.Padding = [6 6 6 6];

            tool.AmpAxes = uiaxes(g);
            tool.AmpAxes.Layout.Row = 1;
            tool.TraceAxes = uiaxes(g);
            tool.TraceAxes.Layout.Row = 2;

            btns = uigridlayout(g, [1 6]);
            btns.Layout.Row = 3;
            btns.Padding = [0 0 0 0];
            btns.ColumnWidth = {110, 110, 115, 120, 160, "1x"};
            uibutton(btns, Text="Lasso select", BackgroundColor=[0.9 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.lassoSelect());
            uibutton(btns, Text="Split bands", BackgroundColor=[0.9 1 0.9], ...
                Tooltip="Detect horizontal amplitude bands (drift-tolerant)", ...
                ButtonPushedFcn=@(~, ~) tool.splitBands());
            uibutton(btns, Text="Snap start -> 0", ...
                ButtonPushedFcn=@(~, ~) tool.snapStart());
            uibutton(btns, Text="Snap stop -> end", ...
                ButtonPushedFcn=@(~, ~) tool.snapStop());
            uibutton(btns, Text="Revert (whole recording)", ...
                ButtonPushedFcn=@(~, ~) tool.revert());
            tool.InfoLabel = uilabel(btns, Text="");

            bandRow = uigridlayout(g, [1 4]);
            bandRow.Layout.Row = 4;
            bandRow.Padding = [0 0 0 0];
            bandRow.ColumnWidth = {200, 160, 200, "1x"};
            tool.BandDropdown = uidropdown(bandRow, Items={'(no bands yet)'}, ...
                Tooltip="Pick a band to see its spikes in the Curation widget", ...
                ValueChangedFcn=@(s, ~) tool.selectBand(s.Value));
            tool.BandMoveBtn = uibutton(bandRow, Text="Band -> new cluster", ...
                Enable="off", ButtonPushedFcn=@(~, ~) tool.moveBand());
            tool.BandMoveAllBtn = uibutton(bandRow, ...
                Text="Out-of-bands -> new clusters", Enable="off", ...
                Tooltip="Give every band except the main one its own cluster", ...
                ButtonPushedFcn=@(~, ~) tool.moveOutOfBands());
            uilabel(bandRow, Text="band 1 = main band", ...
                FontColor=[0.45 0.45 0.45]);
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
            % Colour by detected band when Split bands has run, else by cluster.
            lab = zeros(numel(allIdx), 1);
            if ~isempty(tool.BandIdx)
                [tf, loc] = ismember(allIdx, tool.BandIdx);
                lab(tf) = tool.BandLabels(loc(tf));
            end
            if any(lab > 0)
                k = max(lab);
                pal = lines(k);
                hs = gobjects(1, k);
                for b = 1:k
                    m = lab == b;
                    hs(b) = scatter(ax, tool.PlotX(m), tool.PlotY(m), 4, ...
                        pal(b, :), "filled", MarkerFaceAlpha=0.3);
                end
                if any(lab == 0)
                    scatter(ax, tool.PlotX(lab == 0), tool.PlotY(lab == 0), ...
                        4, [0.7 0.7 0.7], "filled", MarkerFaceAlpha=0.3);
                end
                legend(ax, hs, "band " + string(1:k), Location="best", ...
                    Box="off", AutoUpdate="off");
            else
                for cid = tool.ClusterIds
                    m = s.clusters(allIdx) == cid;
                    scatter(ax, tool.PlotX(m), tool.PlotY(m), 4, ...
                        tool.App.clusterColor(cid), "filled", MarkerFaceAlpha=0.3);
                end
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

        function splitBands(tool)
            %SPLITBANDS Detect horizontal amplitude bands and number them by
            %   size (band 1 = main band). Two passes: first against a static
            %   reference (the overall median), then against a drift trend
            %   computed from the main band only - a running median over ALL
            %   spikes jumps between bands when two are comparable in size,
            %   which smears the ratio density and hides the bands.
            % Only the spikes inside the current [start, stop] window take
            % part - place the bars around an epoch to band-split just that.
            w = tool.App.timeWindowSec();
            inWin = tool.PlotX >= w(1) & tool.PlotX <= w(2);
            t = tool.PlotX(inWin);
            a = tool.PlotY(inWin);
            n = numel(a);
            if n < 50
                tool.InfoLabel.Text = "Too few spikes in the window to find bands.";
                return
            end
            [lab, counts] = tool.bandsFromRatio( ...
                log2(max(a, eps) / max(median(a), eps)));
            trend = tool.bandTrend(t, a, lab == 1, w);
            [lab2, counts2] = tool.bandsFromRatio( ...
                log2(max(a, eps) ./ max(trend, eps)));
            if numel(counts2) >= numel(counts)   % keep the better separation
                lab = lab2;
                counts = counts2;
            end
            tool.BandIdx = tool.PlotIdx(inWin);
            tool.BandLabels = lab;
            k = numel(counts);
            tool.BandDropdown.Items = ["(pick a band)"; ...
                "band " + string(1:k)' + "  (" + string(counts(:)) + " spikes)"];
            tool.BandDropdown.ItemsData = [{[]}; num2cell(1:k)'];
            en = matlab.lang.OnOffSwitchState(k >= 2);
            tool.BandMoveBtn.Enable = en;
            tool.BandMoveAllBtn.Enable = en;
            tool.redrawAmp();
            epoch = "";
            if w(1) > 0 || w(2) < tool.duration()
                epoch = sprintf(" in %.1f-%.1f s", w(1), w(2));
            end
            if k < 2
                tool.InfoLabel.Text = "Only one band found" + epoch + ...
                    " - nothing to split.";
            else
                tool.InfoLabel.Text = sprintf("%d bands found%s (band 1 = " + ...
                    "main). Pick one to inspect, or move the out-of-bands.", ...
                    k, epoch);
            end
        end

        function trend = bandTrend(~, t, a, mask, w)
            %BANDTREND Smoothed running median of a(mask) over time (the
            %   drift of one band), interpolated onto every t. Falls back to
            %   the static median when too few bins have data.
            nb = max(10, min(60, floor(sum(mask) / 100)));
            edges = linspace(w(1), w(2), nb + 1);
            ctr = (edges(1:end-1) + edges(2:end)) / 2;
            med = nan(1, nb);
            for i = 1:nb
                bi = mask & t >= edges(i) & t < edges(i + 1);
                if sum(bi) >= 5
                    med(i) = median(a(bi));
                end
            end
            ok = find(~isnan(med));
            if numel(ok) >= 2
                med = movmean(med, 3, "omitnan");
                tc = min(max(t, ctr(ok(1))), ctr(ok(end)));
                trend = interp1(ctr(ok), med(ok), tc, "linear");
            else
                trend = repmat(median(a(mask)), size(a));
            end
        end

        function [lab, counts] = bandsFromRatio(~, lr)
            %BANDSFROMRATIO Band label per spike from log-amplitude ratios.
            %   Bands = modes of the ratio density, merged until every valley
            %   is deep (density below half the smaller neighbouring peak, a
            %   scale-free test so sparse bands survive next to dense ones),
            %   tiny slivers folded in, then numbered by size (1 = biggest).
            n = numel(lr);
            [f, xi] = ksdensity(lr, NumPoints=512);
            pk = find(islocalmax(f));
            vl = find(islocalmin(f));
            while numel(pk) > 1
                sep = zeros(1, numel(pk) - 1);
                for j = 1:numel(pk) - 1
                    vs = vl(vl > pk(j) & vl < pk(j + 1));
                    sep(j) = min(f(vs)) / min(f(pk(j)), f(pk(j + 1)));
                end
                [worst, wj] = max(sep);
                if worst <= 0.5
                    break
                end
                if f(pk(wj)) < f(pk(wj + 1))
                    pk(wj) = [];
                else
                    pk(wj + 1) = [];
                end
            end
            cuts = zeros(1, numel(pk) - 1);
            for j = 1:numel(pk) - 1
                vs = vl(vl > pk(j) & vl < pk(j + 1));
                [~, vi] = min(f(vs));
                cuts(j) = xi(vs(vi));
            end
            raw = discretize(lr, [-inf, cuts, inf]);
            % Fold tiny slivers into their bigger neighbour.
            minCount = max(20, round(0.002 * n));
            cnt = accumarray(raw(:), 1, [numel(cuts) + 1, 1]);
            while ~isempty(cuts) && min(cnt) < minCount
                [~, i] = min(cnt);
                if i == 1
                    drop = 1;
                elseif i == numel(cnt)
                    drop = numel(cuts);
                elseif cnt(i - 1) >= cnt(i + 1)
                    drop = i - 1;
                else
                    drop = i;
                end
                cuts(drop) = [];
                raw = discretize(lr, [-inf, cuts, inf]);
                cnt = accumarray(raw(:), 1, [numel(cuts) + 1, 1]);
            end
            [counts, gs] = groupcounts(raw(:));
            [counts, order] = sort(counts, "descend");
            gs = gs(order);
            lab = zeros(n, 1);
            for i = 1:numel(gs)
                lab(raw == gs(i)) = i;
            end
        end

        function idx = bandSpikes(tool, b)
            idx = tool.BandIdx(tool.BandLabels == b);
        end

        function selectBand(tool, b)
            %SELECTBAND Show a band's spikes as the Curation selection.
            if isempty(tool.BandLabels) || ~isvalid(tool.Curate)
                return
            end
            if isempty(b)
                tool.Curate.selectGlobal([]);   % (pick a band) clears it
                return
            end
            idx = tool.bandSpikes(b);
            tool.Curate.selectGlobal(idx);
            tool.InfoLabel.Text = sprintf("Band %d: %d spikes " + ...
                "(shown in the Curation widget).", b, numel(idx));
        end

        function moveBand(tool)
            %MOVEBAND Split the chosen band off into a new cluster.
            b = tool.BandDropdown.Value;
            if isempty(tool.BandLabels) || isempty(b)
                tool.InfoLabel.Text = "Pick a band in the dropdown first.";
                return
            end
            idx = tool.bandSpikes(b);
            newId = max(tool.App.Spikes.clusterIds) + 1;
            tool.App.applySplit(idx(:), newId);
            tool.adoptIntoCuration(newId);
            tool.clearBands(sprintf("Band %d -> new cluster %d (%d spikes).", ...
                b, newId, numel(idx)));
        end

        function moveOutOfBands(tool)
            %MOVEOUTOFBANDS Every band except the main one -> its own cluster.
            if isempty(tool.BandLabels) || max(tool.BandLabels) < 2
                tool.InfoLabel.Text = "No out-of-bands (run Split bands first).";
                return
            end
            newIds = [];
            moved = strings(0, 1);
            for b = 2:max(tool.BandLabels)
                idx = tool.bandSpikes(b);
                if isempty(idx)
                    continue
                end
                newId = max(tool.App.Spikes.clusterIds) + 1;
                tool.App.applySplit(idx(:), newId);
                newIds(end + 1) = newId; %#ok<AGROW>
                moved(end + 1) = sprintf("band %d -> cluster %d", b, newId); %#ok<AGROW>
            end
            tool.adoptIntoCuration(newIds);
            tool.clearBands("Moved " + join(moved, ", ") + ".");
        end

        function adoptIntoCuration(tool, newIds)
            %ADOPTINTOCURATION Keep the split-off clusters visible in the
            %   Curation widget (and therefore in this amplitude plot).
            if isempty(newIds)
                return
            end
            if ~isempty(tool.Curate) && isvalid(tool.Curate)
                tool.Curate.adoptCluster(newIds);
                tool.ClusterIds = tool.Curate.shownClusters();
            else
                tool.ClusterIds = unique([tool.ClusterIds, newIds]);
            end
        end

        function clearBands(tool, msg)
            %CLEARBANDS Band labels are stale after any cluster edit.
            tool.BandIdx = [];
            tool.BandLabels = [];
            tool.BandDropdown.Items = {'(no bands yet)'};
            tool.BandDropdown.ItemsData = [];
            tool.BandMoveBtn.Enable = "off";
            tool.BandMoveAllBtn.Enable = "off";
            tool.redrawAmp();
            tool.InfoLabel.Text = string(msg);
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
