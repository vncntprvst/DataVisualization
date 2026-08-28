classdef RecoverSpikesTool < handle
    %RECOVERSPIKESTOOL Re-detect spikes the sorter missed, by template match.
    %   RecoverSpikesTool(app, curate) opens a companion window (Curation ->
    %   Tools) that finds events on the raw trace matching a cluster's mean
    %   waveform but not already in the sorting, so you can recover e.g. the
    %   low-amplitude spikes a unit had early in a recording.
    %
    %   For the target cluster it slides the mean waveform over the data channel
    %   within a time window, keeps threshold-crossing events whose shape
    %   correlates with the template above a cutoff and that are not within the
    %   refractory period of any existing spike, previews them on the amplitude
    %   and raw-trace panels, and (on Accept) adds them to a cluster.
    %
    %   Raw trace: scroll to pan, Ctrl+scroll (or up/down arrows) to zoom. Drag
    %   the red threshold line to set the detection level (the "thr uV" field
    %   tracks it, and editing the field moves the line). Accept adds the
    %   candidates to the target cluster, or to a brand-new cluster when the
    %   "new cluster" box is ticked.
    %
    %   See also CurationTool, SpikeVisualizationApp.

    properties (Access = private)
        App
        Curate

        UIFigure               matlab.ui.Figure
        AmpAxes                matlab.ui.control.UIAxes
        TraceAxes              matlab.ui.control.UIAxes
        TargetDropdown         matlab.ui.control.DropDown
        StartField             matlab.ui.control.NumericEditField
        StopField              matlab.ui.control.NumericEditField
        ThreshField            matlab.ui.control.NumericEditField
        CorrField              matlab.ui.control.NumericEditField
        DedupField             matlab.ui.control.NumericEditField
        NewClusterBox          matlab.ui.control.CheckBox
        InfoLabel              matlab.ui.control.Label

        ThresholdLine = images.roi.Line.empty   % draggable detection level
        TraceSlider matlab.ui.control.Slider    % scrub the trace excerpt
        SuppressLineEvents logical = false
        TraceCenterSec double = NaN              % raw-trace pan/zoom state
        TraceHalfSec double = 0.125

        CandTimes double = []   % candidate spike sample times
        CandAmp double = []     % candidate peak-to-peak, microvolts
    end

    properties (Constant, Access = private)
        RefractoryMs = 1
    end

    methods
        function tool = RecoverSpikesTool(app, curate)
            arguments
                app (1,1) SpikeVisualizationApp
                curate (1,1) CurationTool
            end
            tool.App = app;
            tool.Curate = curate;
            tool.buildUI();
            tool.refreshPlots();
            tool.setInfo("Set the target cluster, window and thresholds, then Detect.");
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

    methods (Access = private)
        function buildUI(tool)
            tool.UIFigure = uifigure(Name="Recover missing spikes", ...
                Position=[180 120 1120 640], ...
                WindowScrollWheelFcn=@(~, e) tool.onScroll(e), ...
                WindowKeyPressFcn=@(~, e) tool.onKey(e), ...
                CloseRequestFcn=@(~, ~) delete(tool));
            g = uigridlayout(tool.UIFigure, [5 1]);
            g.RowHeight = {"1x", "1x", 22, 32, 18};
            g.RowSpacing = 6;
            g.Padding = [6 6 6 6];

            tool.AmpAxes = uiaxes(g);
            tool.AmpAxes.Layout.Row = 1;
            tool.TraceAxes = uiaxes(g);
            tool.TraceAxes.Layout.Row = 2;
            tool.TraceSlider = uislider(g, Limits=[0 1], Value=0, ...
                MajorTicks=[], MinorTicks=[], ...
                ValueChangedFcn=@(~, e) tool.onTraceSlider(e));
            tool.TraceSlider.Layout.Row = 3;

            ctl = uigridlayout(g, [1 16]);
            ctl.Layout.Row = 4;
            ctl.Padding = [0 0 0 0];
            ctl.ColumnSpacing = 5;
            ctl.ColumnWidth = {"fit", 55, "fit", 55, "fit", 55, "fit", 55, ...
                "fit", 55, "fit", 55, 65, 65, 100, 55};

            ids = tool.Curate.shownClusters();
            uilabel(ctl, Text="cluster", HorizontalAlignment="right");
            tool.TargetDropdown = uidropdown(ctl, Items=string(ids(:)), ...
                ItemsData=ids(:), Value=tool.defaultTarget(ids), ...
                ValueChangedFcn=@(~, ~) tool.onTargetChanged());
            w = tool.App.timeWindowSec();
            uilabel(ctl, Text="start s", HorizontalAlignment="right");
            tool.StartField = uieditfield(ctl, "numeric", Value=w(1));
            uilabel(ctl, Text="stop s", HorizontalAlignment="right");
            tool.StopField = uieditfield(ctl, "numeric", Value=w(2));
            uilabel(ctl, Text="thr uV", HorizontalAlignment="right");
            tool.ThreshField = uieditfield(ctl, "numeric", ...
                Value=tool.defaultThreshold(tool.defaultTarget(ids)), ...
                Tooltip="Signed detection level: negative catches downward " + ...
                "deflections, positive upward (the line sits where it detects)", ...
                ValueChangedFcn=@(~, ~) tool.onThreshFieldChanged());
            uilabel(ctl, Text="min corr", HorizontalAlignment="right");
            tool.CorrField = uieditfield(ctl, "numeric", Value=0.8, ...
                Limits=[0 1]);
            uilabel(ctl, Text="dedup ms", HorizontalAlignment="right");
            tool.DedupField = uieditfield(ctl, "numeric", ...
                Value=tool.RefractoryMs, Limits=[0 Inf], ...
                Tooltip="Drop candidates within this window of any existing " + ...
                "spike (in any cluster), to avoid duplicates");
            uibutton(ctl, Text="Detect", BackgroundColor=[0.9 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.detect());
            uibutton(ctl, Text="Accept", BackgroundColor=[0.85 1 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.accept());
            tool.NewClusterBox = uicheckbox(ctl, Text="new cluster");
            uibutton(ctl, Text="Clear", ...
                ButtonPushedFcn=@(~, ~) tool.clearCandidates());

            tool.InfoLabel = uilabel(g, Text="");
            tool.InfoLabel.Layout.Row = 5;
        end

        function onTraceSlider(tool, event)
            %ONTRACESLIDER Scrub the trace excerpt through the recording.
            s = tool.App.Spikes;
            tool.TraceCenterSec = event.Value * s.numSamples / s.samplingRate;
            tool.plotTrace();
        end

        function cid = defaultTarget(tool, ids)
            counts = arrayfun(@(c) sum(tool.App.Spikes.clusters == c), ids);
            [~, k] = max(counts);
            cid = ids(k);
        end

        function thr = defaultThreshold(tool, cid)
            % ~40% of the template's dominant peak. SIGNED: the sign picks
            % the detection side (negative = downward deflections).
            idx = tool.App.Spikes.clusters == cid;
            uv = tool.uvScale();
            template = mean(double(tool.App.Spikes.waveforms(idx, :)), 1) * uv;
            [m, tPeak] = max(abs(template));
            thr = round(0.4 * m) * sign(template(tPeak));
            if ~isfinite(thr) || thr == 0
                thr = -50;
            end
        end

        function cid = targetCluster(tool)
            cid = tool.TargetDropdown.Value;
        end

        function s = uvScale(tool)
            s = 1;
            if isfield(tool.App.Spikes, "uvPerADC")
                s = tool.App.Spikes.uvPerADC;
            end
        end

        function onTargetChanged(tool)
            tool.ThreshField.Value = tool.defaultThreshold(tool.targetCluster());
            tool.TraceCenterSec = tool.StartField.Value;
            tool.clearCandidates();
        end

        % ------------------------------------------------- trace pan / zoom
        function d = traceDurSec(tool)
            s = tool.App.Spikes;
            d = s.numSamples / s.samplingRate;
        end

        function onScroll(tool, event)
            if isnan(tool.TraceCenterSec)
                return
            end
            if any(strcmp(tool.UIFigure.CurrentModifier, "control"))
                tool.zoomTrace(1.25 ^ event.VerticalScrollCount);
            else
                tool.panTrace(0.2 * event.VerticalScrollCount);
            end
        end

        function onKey(tool, event)
            if isnan(tool.TraceCenterSec)
                return
            end
            switch event.Key
                case "leftarrow"
                    tool.panTrace(-0.25);
                case "rightarrow"
                    tool.panTrace(0.25);
                case {"uparrow", "add", "equal"}
                    tool.zoomTrace(1 / 1.25);
                case {"downarrow", "subtract", "hyphen"}
                    tool.zoomTrace(1.25);
            end
        end

        function zoomTrace(tool, factor)
            tool.TraceHalfSec = min(5, max(0.005, tool.TraceHalfSec * factor));
            tool.clampCenter();
            tool.plotTrace();
        end

        function panTrace(tool, fractionOfWindow)
            tool.TraceCenterSec = tool.TraceCenterSec + ...
                fractionOfWindow * 2 * tool.TraceHalfSec;
            tool.clampCenter();
            tool.plotTrace();
        end

        function clampCenter(tool)
            dur = tool.traceDurSec();
            h = tool.TraceHalfSec;
            tool.TraceCenterSec = min(max(tool.TraceCenterSec, h), max(h, dur - h));
        end

        % ------------------------------------------------- threshold line
        function drawThresholdLine(tool)
            ax = tool.TraceAxes;
            tool.SuppressLineEvents = true;
            if ~isempty(tool.ThresholdLine) && isvalid(tool.ThresholdLine)
                delete(tool.ThresholdLine);
            end
            tool.ThresholdLine = images.roi.Line.empty;
            if isvalid(ax)
                yv = tool.ThreshField.Value;   % signed: drawn where it detects
                xl = xlim(ax);
                roi = images.roi.Line(ax, Position=[xl(1) yv; xl(2) yv], ...
                    Color=[0.85 0 0], LineWidth=2.5, MarkerSize=0.1, ...
                    Deletable=false);
                addlistener(roi, "MovingROI", @(src, ~) tool.constrainThreshLine(src));
                addlistener(roi, "ROIMoved", @(~, ~) tool.onThreshLineMoved());
                tool.ThresholdLine = roi;
            end
            tool.SuppressLineEvents = false;
        end

        function constrainThreshLine(tool, roi)
            if tool.SuppressLineEvents || ~isvalid(roi)
                return
            end
            ym = mean(roi.Position(:, 2));
            xl = xlim(tool.TraceAxes);
            roi.Position = [xl(1) ym; xl(2) ym];   % keep it horizontal
        end

        function onThreshLineMoved(tool)
            if tool.SuppressLineEvents || isempty(tool.ThresholdLine) ...
                    || ~isvalid(tool.ThresholdLine)
                return
            end
            ym = mean(tool.ThresholdLine.Position(:, 2));
            thr = round(ym);
            if thr == 0   % keep a definite side; 0 has no polarity
                if ym < 0
                    thr = -1;
                else
                    thr = 1;
                end
            end
            tool.ThreshField.Value = thr;
            SpikeVisualizationApp.logEvent(sprintf( ...
                "recover thr drag: ym=%.1f -> thr=%d", ym, thr));
            tool.detect();
        end

        function onThreshFieldChanged(tool)
            tool.drawThresholdLine();
            tool.detect();
        end

        function detect(tool)
            s = tool.App.Spikes;
            cid = tool.targetCluster();
            idx = s.clusters == cid;
            if nnz(idx) < 5
                tool.setInfo("Target cluster has too few spikes for a template.");
                return
            end
            uvScl = tool.uvScale();
            pre = s.waveformWindow(1);
            post = s.waveformWindow(2);
            template = mean(double(s.waveforms(idx, :)), 1) * uvScl;   % microvolts
            % The threshold's SIGN picks the detection side (drag the line
            % across zero to hunt the opposite lobe); align candidates on
            % the template's extremum on that side.
            detPol = sign(tool.ThreshField.Value);
            if detPol == 0
                [~, tp] = max(abs(template));
                detPol = sign(template(tp));
                if detPol == 0
                    detPol = -1;
                end
            end
            [~, tPeak] = max(template * detPol);
            offset = tPeak - (pre + 1);

            fs = s.samplingRate;
            first = max(1, round(tool.StartField.Value * fs));
            last = min(s.numSamples, round(tool.StopField.Value * fs));
            if last - first < pre + post + 2
                tool.setInfo("Window too short.");
                return
            end
            tool.setInfo("Detecting...");
            drawnow;
            uv = double(tool.App.traceSamples(first, last));   % microvolts
            if isempty(uv)
                tool.setInfo("Raw trace unavailable.");
                return
            end

            refr = round(fs * tool.RefractoryMs / 1000);
            [~, locs] = findpeaks(uv * detPol, ...
                MinPeakHeight=abs(tool.ThreshField.Value), MinPeakDistance=refr);
            if numel(locs) > 200000
                tool.setInfo(sprintf("%d raw crossings at thr %d uV - too " + ...
                    "many to score. Move the threshold further from zero.", ...
                    numel(locs), tool.ThreshField.Value));
                return
            end
            spkTimes = first - 1 + locs - offset;   % candidate spike sample times

            % Keep candidates whose snippet matches the template shape.
            score = nan(numel(spkTimes), 1);
            amp = nan(numel(spkTimes), 1);
            for k = 1:numel(spkTimes)
                lo = spkTimes(k) - first + 1 - pre;
                hi = spkTimes(k) - first + 1 + post;
                if lo >= 1 && hi <= numel(uv)
                    snip = uv(lo:hi);
                    score(k) = corr(snip(:), template(:));
                    amp(k) = max(snip) - min(snip);
                end
            end
            keep = score >= tool.CorrField.Value;
            spkTimes = spkTimes(keep);
            amp = amp(keep);

            % Drop anything already sorted (within the dedup window of any spike
            % in any cluster), so accepting never creates duplicates.
            dedup = max(1, round(fs * tool.DedupField.Value / 1000));
            existing = unique(s.spikeTimes);
            if numel(existing) >= 2 && ~isempty(spkTimes)
                nearest = interp1(existing, existing, spkTimes, "nearest", "extrap");
                far = abs(nearest - spkTimes) >= dedup;
                spkTimes = spkTimes(far);
                amp = amp(far);
            end

            tool.CandTimes = spkTimes;
            tool.CandAmp = amp;
            if ~isempty(spkTimes)
                tool.TraceCenterSec = spkTimes(1) / fs;   % focus on the recovered spikes
            end
            tool.refreshPlots();
            tool.setInfo(sprintf("%d new candidates in %.1f-%.1f s " + ...
                "(thr %.0f uV, corr>=%.2f). Review, then Accept.", numel(spkTimes), ...
                tool.StartField.Value, tool.StopField.Value, ...
                tool.ThreshField.Value, tool.CorrField.Value));
        end

        function accept(tool)
            if isempty(tool.CandTimes)
                tool.setInfo("Nothing to accept - Detect first.");
                return
            end
            n = numel(tool.CandTimes);
            if tool.NewClusterBox.Value
                cid = max(tool.App.Spikes.clusterIds) + 1;   % fresh cluster id
                suffix = " (new cluster)";
            else
                cid = tool.targetCluster();
                suffix = "";
            end
            tool.App.addSpikes(tool.CandTimes, cid);
            tool.CandTimes = [];
            tool.CandAmp = [];
            tool.refreshDropdown();
            tool.refreshPlots();
            tool.setInfo(sprintf("Accepted %d spikes into cluster %d%s.", ...
                n, cid, suffix));
        end

        function refreshDropdown(tool)
            % Keep the target list in step with the sorting (e.g. after a new
            % cluster was added), preserving the current selection if possible.
            ids = unique(tool.App.Spikes.clusters);
            keep = tool.TargetDropdown.Value;
            tool.TargetDropdown.Items = string(ids(:));
            tool.TargetDropdown.ItemsData = ids(:);
            if ismember(keep, ids)
                tool.TargetDropdown.Value = keep;
            end
        end

        function clearCandidates(tool)
            tool.CandTimes = [];
            tool.CandAmp = [];
            tool.refreshPlots();
        end

        function refreshPlots(tool)
            tool.plotAmplitude();
            tool.plotTrace();
        end

        function plotAmplitude(tool)
            ax = tool.AmpAxes;
            cla(ax);
            hold(ax, "on");
            s = tool.App.Spikes;
            cid = tool.targetCluster();
            idx = s.clusters == cid;
            scatter(ax, s.spikeTimes(idx) / s.samplingRate, ...
                s.amplitudePP(idx) * tool.uvScale(), 5, ...
                tool.App.clusterColor(cid), "filled", MarkerFaceAlpha=0.3);
            if ~isempty(tool.CandTimes)
                scatter(ax, tool.CandTimes / s.samplingRate, tool.CandAmp, 12, ...
                    [0 0 0], "x", LineWidth=0.75);
            end
            hold(ax, "off");
            title(ax, sprintf("Amplitude vs time  (cluster %d + %d candidates)", ...
                cid, numel(tool.CandTimes)));
            xlabel(ax, "Time (s)");
            ylabel(ax, "Peak-to-peak (\muV)");
            xlim(ax, [0 s.numSamples / s.samplingRate]);
        end

        function plotTrace(tool)
            ax = tool.TraceAxes;
            cla(ax);
            s = tool.App.Spikes;
            fs = s.samplingRate;
            if isnan(tool.TraceCenterSec)
                if ~isempty(tool.CandTimes)
                    tool.TraceCenterSec = tool.CandTimes(1) / fs;
                else
                    tool.TraceCenterSec = tool.StartField.Value;
                end
            end
            [tSec, uv] = tool.App.traceExcerpt(tool.TraceCenterSec, tool.TraceHalfSec);
            if isempty(tSec)
                tool.ThresholdLine = images.roi.Line.empty;
                title(ax, "Raw trace unavailable");
                return
            end
            plot(ax, tSec, uv, Color=[0 0 0 0.8], LineWidth=0.1);
            hold(ax, "on");
            yb = double(min(uv));
            cid = tool.targetCluster();
            ex = s.spikeTimes(s.clusters == cid) / fs;
            ex = ex(ex >= tSec(1) & ex <= tSec(end));
            plot(ax, ex, repmat(yb, size(ex)), "^", Color=tool.App.clusterColor(cid), ...
                MarkerFaceColor=tool.App.clusterColor(cid), MarkerSize=5, ...
                LineStyle="none");
            if ~isempty(tool.CandTimes)
                cd = tool.CandTimes / fs;
                cd = cd(cd >= tSec(1) & cd <= tSec(end));
                plot(ax, cd, repmat(yb, size(cd)), "x", Color=[0 0 0], ...
                    MarkerSize=7, LineWidth=1, LineStyle="none");
            end
            hold(ax, "off");
            title(ax, "Raw trace (^ existing, x recovered)  -  " + ...
                "scroll to pan, Ctrl+scroll to zoom, drag the red line for thr");
            xlabel(ax, "Time (s)");
            ylabel(ax, "Voltage (\muV)");
            xlim(ax, [tSec(1) tSec(end)]);
            tool.drawThresholdLine();
            if ~isempty(tool.TraceSlider) && isvalid(tool.TraceSlider)
                dur = s.numSamples / fs;
                tool.TraceSlider.Value = ...
                    min(1, max(0, tool.TraceCenterSec / dur));
            end
        end

        function setInfo(tool, msg)
            if ~isempty(tool.InfoLabel) && isvalid(tool.InfoLabel)
                tool.InfoLabel.Text = string(msg);
            end
        end
    end
end
