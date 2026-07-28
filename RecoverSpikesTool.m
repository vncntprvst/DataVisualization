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
    %   and raw-trace panels, and (on Accept) adds them to the cluster.
    %
    %   Detection threshold (microvolts) and correlation cutoff are adjustable;
    %   lower the threshold to catch smaller spikes, raise the correlation to
    %   stay specific to the cluster's shape.
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
        InfoLabel              matlab.ui.control.Label

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
    end

    methods (Access = private)
        function buildUI(tool)
            tool.UIFigure = uifigure(Name="Recover missing spikes", ...
                Position=[180 120 1000 640]);
            g = uigridlayout(tool.UIFigure, [4 1]);
            g.RowHeight = {"1x", "1x", 34, 18};
            g.Padding = [6 6 6 6];

            tool.AmpAxes = uiaxes(g);
            tool.AmpAxes.Layout.Row = 1;
            tool.TraceAxes = uiaxes(g);
            tool.TraceAxes.Layout.Row = 2;

            ctl = uigridlayout(g, [1 12]);
            ctl.Layout.Row = 3;
            ctl.Padding = [0 0 0 0];
            ctl.ColumnWidth = {"fit", 60, "fit", 55, "fit", 55, "fit", 55, ...
                "fit", 90, 80, 80};

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
                Value=tool.defaultThreshold(tool.defaultTarget(ids)));
            uilabel(ctl, Text="min corr", HorizontalAlignment="right");
            tool.CorrField = uieditfield(ctl, "numeric", Value=0.8, ...
                Limits=[0 1]);
            uibutton(ctl, Text="Detect", BackgroundColor=[0.9 0.9 1], ...
                ButtonPushedFcn=@(~, ~) tool.detect());
            uibutton(ctl, Text="Accept", BackgroundColor=[0.85 1 0.85], ...
                ButtonPushedFcn=@(~, ~) tool.accept());
            uibutton(ctl, Text="Clear", ...
                ButtonPushedFcn=@(~, ~) tool.clearCandidates());

            tool.InfoLabel = uilabel(g, Text="");
            tool.InfoLabel.Layout.Row = 4;
        end

        function cid = defaultTarget(tool, ids)
            counts = arrayfun(@(c) sum(tool.App.Spikes.clusters == c), ids);
            [~, k] = max(counts);
            cid = ids(k);
        end

        function thr = defaultThreshold(tool, cid)
            % ~40% of the template's dominant peak, a starting point to tune.
            idx = tool.App.Spikes.clusters == cid;
            uv = tool.uvScale();
            template = mean(double(tool.App.Spikes.waveforms(idx, :)), 1) * uv;
            thr = round(0.4 * max(abs(template)));
            if ~isfinite(thr) || thr <= 0
                thr = 50;
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
            tool.clearCandidates();
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
            [~, tPeak] = max(abs(template));
            polarity = sign(template(tPeak));
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
            [~, locs] = findpeaks(uv * polarity, ...
                MinPeakHeight=tool.ThreshField.Value, MinPeakDistance=refr);
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

            % Drop anything already sorted (within the refractory of any spike).
            existing = unique(s.spikeTimes);
            if numel(existing) >= 2 && ~isempty(spkTimes)
                nearest = interp1(existing, existing, spkTimes, "nearest", "extrap");
                far = abs(nearest - spkTimes) >= refr;
                spkTimes = spkTimes(far);
                amp = amp(far);
            end

            tool.CandTimes = spkTimes;
            tool.CandAmp = amp;
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
            cid = tool.targetCluster();
            tool.App.addSpikes(tool.CandTimes, cid);
            tool.CandTimes = [];
            tool.CandAmp = [];
            tool.refreshPlots();
            tool.setInfo(sprintf("Accepted %d spikes into cluster %d.", n, cid));
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
            % Focus on the first candidate (else the window start) so recovered
            % spikes are visible.
            if ~isempty(tool.CandTimes)
                centre = tool.CandTimes(1) / fs;
            else
                centre = tool.StartField.Value;
            end
            [tSec, uv] = tool.App.traceExcerpt(centre, 0.25);
            if isempty(tSec)
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
            title(ax, "Raw trace (^ existing, x recovered)");
            xlabel(ax, "Time (s)");
            ylabel(ax, "Voltage (\muV)");
            axis(ax, "tight");
        end

        function setInfo(tool, msg)
            if ~isempty(tool.InfoLabel) && isvalid(tool.InfoLabel)
                tool.InfoLabel.Text = string(msg);
            end
        end
    end
end
