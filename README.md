# SpikeVisualizationApp

Interactive viewer and curation tool for sorted extracellular spikes, focused on
single-electrode recordings exported to the **Phy / KiloSort** format. It is the
modern `uifigure` replacement for the legacy GUIDE `SpikeVisualizationGUI`
(retired; recoverable from git history).

Requires **MATLAB R2025a** with the Statistics & Machine Learning, Signal
Processing and Image Processing toolboxes, plus
[`npy-matlab`](https://github.com/kwikteam/npy-matlab) (`readNPY`/`writeNPY`) on
the path.

## Files

| File | Purpose |
|------|---------|
| `SpikeVisualizationApp.m` | Main app: cluster list, waveforms, mean±SD, amplitude drift, ISI (log), ACG, PC features, and a scrollable raw-trace view with spike markers, zoom, and threshold lines. Cluster labelling, merging and saving. |
| `CurationTool.m` | Companion window with three linked panels (ISI, PC features, waveforms). Select spikes by circling ISI bars, lassoing PC features, or drawing a line across the waveforms, then split or merge. |
| `TimeWindowTool.m` | Companion (Curation → Tools) that restricts curation to a draggable time band — for splitting units that swap amplitude over time on one electrode. |
| `RecoverSpikesTool.m` | Companion (Curation → Tools) that re-detects spikes the sorter missed by template-matching a cluster's mean waveform against the raw trace. |
| `PETHTool.m` | Companion (main app → Tools) that shows peri-event time histograms of the selected clusters aligned to loaded behavioural events. |
| `LoadEvents.m` | Reader for an events file (`events.mat`, or long `events.tsv`/`.csv`) of alignment times in recording seconds. |
| `DatasetBrowser.m` | Companion (main app → Tools, or standalone) that browses a dataset table by unit or recording, filters/sorts it, and opens an entry's sorting (and unit) in the app. |
| `LoadDataset.m` | Reader for a dataset table (`dataset.tsv`/`.csv`, one row per unit) — returns the table plus column info. |
| `LoadSpikesPhy.m` | Reader for a Phy export folder (`params.py`, `spike_*.npy`, `pc_features.npy`, `cluster_*.tsv`), extracting one waveform per spike from the raw data channel. |
| `LoadSpikesOnline.m` | Reader for an online-sorted REX session (`<tag>_REX.mat`) — spike times only, no raw trace or waveforms. |
| `addons/*.m` | Drop-in scripts exposed under Options → Run script. |
| `addons/add_cs_cluster.m` | Injects the P-sort complex-spike train into a Phy folder as its own cluster (see below). |

## Complex spikes (P-sort)

Complex spikes are a separate train and are not part of the ordinary sorting. To
see them, inject the P-sort complex-spike output (`purkinje/purkinje_*.npz`,
found automatically in the session tree) into the Phy folder as a **new, separate
cluster** — it is never pooled with a simple-spike cluster:

```matlab
add_cs_cluster("D:\...\S113L4A5_12540\tdc\phy_postmerge")
```

After that, the app (and Phy) load the complex spikes as a normal extra cluster
with its own waveforms, ISI and markers. The script backs up the originals as
`*.preCS.npy`, is idempotent (writes a `.cs_injected.json` marker), and can be
undone by restoring the backups. There is no complex-spike-specific code in the
app itself.

## Opening a sorting

The sortings you open in Phy with, e.g.

```powershell
.\open_phy.ps1 S113L4A5_12540 tdc post
```

open the folder `...\curation_test\S113L4A5_12540\tdc\phy_postmerge`. Open the
same folder here:

```matlab
SpikeVisualizationApp("D:\CB paper\1-analyses\spikesort_cluster\curation_test\S113L4A5_12540\tdc\phy_postmerge")
```

Or call `SpikeVisualizationApp` with no argument and pick the folder in the
dialog. These exports store `data.dat` as 2 channels where **ch0 is the real
signal** and ch1 is a synthetic duplicate; the loader reads ch0 by default
(`LoadSpikesPhy(dir, DataChannel=0)`).

### Units

Voltages are shown in **microvolts** when the export declares a gain. For these
exports `gain_to_uV` is the number of ADC counts per microvolt, so the app uses
`µV = ADC / gain_to_uV` (≈334 µV median peak-to-peak here, ≈18 µV RMS noise). If
your rig uses the opposite convention, invert `uvPerADC` in `normalizeModel`.

## The panels

The **Clusters** selector on the left is a table with a **colour swatch** per
cluster (matching its colour in the plots), its id, spike count, and label;
select one or more rows.

The plot panels are Mean ± SD (top-left) with the Waveforms overlay beneath it
(its title reports how many of the total are drawn), the ISI histogram (log x
axis) top-middle, Amplitude-vs-time (drift) top-right, a per-cluster PCA feature
scatter, and a **correlogram matrix** — all for the selected cluster(s) — plus a
raw-data excerpt. The cluster table can be sorted by ID, spike count, or label
(the **sort** dropdown in its header).

The correlogram matrix shows one small plot per ordered pair of selected
clusters: the diagonal is each cluster's autocorrelogram (in its colour) and the
off-diagonals are the cross-correlograms (row = reference at lag 0, so the
transposed cell is the reference-flipped view). Up to 4 clusters are shown.

## Edit menu, Reload, last folder

* **Edit → Undo / Redo** (Ctrl+Z / Ctrl+Y) step through curation actions
  (split, merge, realign, relabel).
* **Reload** re-reads the current Phy folder from disk (handy after running a
  script that changed it).
* The **Load Phy folder** dialog reopens at the last folder you loaded.

## Options → Run script (add-ons)

Drop-in scripts placed in the `addons/` folder appear under **Options → Run
script** by a readable name (a `% AddonName: ...` tag in the file, else the
prettified file name). Selecting one runs it on the current Phy folder and then
reloads. `add_cs_cluster` (Add CS cluster) ships as the first such add-on.

### Raw-data excerpt

* **Zoom / pan the time axis:** `Ctrl` + mouse wheel zooms; the wheel alone pans.
  Left/Right arrows pan, Up/Down arrows zoom. The axes zoom/pan toolbar works too
  (it reloads the excerpt for the new limits). The slider jumps anywhere in the
  recording, and **Prev/Next spike** jump the window (same length) to the
  previous/next spike of the selected cluster(s). The voltage axis stays fixed
  while you pan.
* **Threshold lines:** *Add threshold* drops a draggable horizontal red line
  (drag it to a voltage level; press *Delete* to remove it, or *Clear
  thresholds*). Every spike whose waveform crosses **all** threshold lines is
  drawn red on the waveform and mean panels and red-outlined in the amplitude
  and PC feature panels — a quick amplitude-window discriminator. The top-of-app
  dropdown chooses whether this highlights spikes **in the current window**
  (default) or **in the whole recording**.

## Curating: split and merge

1. Select one or more clusters and click **Curation...**. The window has three
   linked panels — ISI histogram, PC features (joint PCA), and waveforms — each
   overlaying every selected cluster in its own colour, so you can compare two
   clusters before merging. The current selection is shown in black. The PC
   projection is chosen with the **PC view** dropdown above the plots, and the
   controls are grouped into **Selection**, **Clustering**, and **Actions**.
2. Select spikes any of these ways (the selection shows in all three panels):
   * **Circle ISI bars** – freehand over the bars close to 0 in the (log) ISI
     histogram; selects the spikes taking part in those short intervals.
   * **Lasso PC features** – freehand in the scatter (temporal PCA of the ch0
     waveform, the same features Phy's FeatureView shows); refines the current
     selection or selects afresh.
   * **Line select (add)** – draw a line across the waveform panel; selects
     **every** shown spike whose waveform crosses that line, not only the ones
     drawn. Click it again to place more lines — a spike must cross **all** of
     them (a corridor).
   * **Find clusters** – automatically split the shown spikes across the PC
     features. It recolours the scatter and waveforms by the found sub-clusters
     (and auto-switches to the PC projection that separates them best), then you
     pick one in the **Found** dropdown to select it. Method (k-means / Gaussian
     mixture) and number of clusters (auto-estimated, or 2–8) are set in the
     widget's **Settings** menu.
   * **Clear lines/found** resets placed lines and any found clustering.
3. **Move selected → [destination]** – reassign exactly the selected spikes to
   the cluster chosen in the dropdown (any existing cluster, or **(new cluster)**
   to split them off). This is how you merge just the selected spikes into
   another cluster, in one step.
4. **Merge shown** – merge **all** the shown clusters into one (optionally plus
   one more, not currently shown, from the "also merge with…" dropdown). Clusters
   merge into the lowest id. This ignores the selection — use *Move selected* to
   act on selected spikes only.

### Splitting units that cross in amplitude (Curation → Tools → Time-window split)

For a single electrode where two units swap amplitude rank over time, open the
**Time-window split** tool. It shows amplitude-vs-time for the curated clusters
with two draggable vertical bars (start/stop) and a raw-trace strip at the bar
you last moved. Moving a bar restricts the Curation widget to the clusters'
spikes **within that time band**, so you can isolate an interval where the two
amplitudes are separated, split there, then step to the next band and stitch.
Buttons snap the start to 0, the stop to the end, or revert to the whole
recording. The tool and the Curation widget are linked live: a selection in
Curation is highlighted on the amplitude plot, and **Lasso select** there sets
the Curation selection — no re-applying or reopening.

### Recovering spikes the sorter missed (Curation → Tools → Recover missing spikes)

To recover spikes a unit lost — e.g. its low-amplitude spikes early in a
recording — open **Recover missing spikes**. It takes the target cluster's mean
waveform as a template, finds threshold-crossing events on the raw trace within a
time window whose shape correlates with the template above a cutoff and that are
**not** within the refractory period of any existing spike, and previews them
(black ×) on the amplitude-vs-time and raw-trace panels. Lower the detection
threshold (µV) to catch smaller spikes; raise the correlation to stay specific to
the cluster's shape. On the raw trace, **scroll to pan, Ctrl+scroll (or up/down
arrows) to zoom**, and **drag the red line** to set the detection threshold (the
`thr uV` field tracks it, and editing the field moves the line). **Accept** adds
them to the target cluster, or — with **new cluster** ticked — to a brand-new
cluster (both undoable). Because the already-sorted spikes are deduped out, in a
window where the unit was under-sorted the recovered events are its missing spikes
— they trace the unit's amplitude drift straight into its sorted spikes.

Recovered candidates are de-duplicated against **all existing spikes in every
cluster** (not just the target): any candidate within the 1 ms refractory of an
already-sorted spike is dropped, so accepting them — into the target or a new
cluster — never creates duplicates.

In the main window, **Merge selected** merges the clusters selected in the list,
**Realign selected** shifts each selected cluster's spike times so its mean
waveform peak lines up at the same sample, then re-extracts the waveforms (use it
when a cluster is a time-shifted duplicate of another), and **Discard selected**
removes the selected cluster(s)' spikes from the sorting entirely (with a confirm;
undoable, and written to disk only when you Save). Cluster labels are **SU / MU /
Noise / Unsorted / Other** (or type any text directly in the Label column).

## Options menu

The **Options** menu (top of the window) has checkable settings:
* *Threshold highlights spikes* — whether the raw-trace threshold highlights
  spikes **in the current window** (default) or **in the whole recording**.
* *Realign to* — the feature **Realign selected** aligns to: the **largest peak**
  (default), the **trough**, or the **first peak** (useful when a later peak is
  larger than the physiologically relevant first one).

## Cluster labels

The *Cluster labels* table lists **every** cluster (including freshly split
ones); Label and Note are editable. Labels are optional — the cluster
*assignment* of every spike is always saved regardless of whether it is labelled.

All curation steps are scriptable for reproducible workflows:

```matlab
app  = SpikeVisualizationApp(phyDir);
tool = CurationTool(app, clusterId);
tool.selectIsiRange([0.2 2]);               % circle ISI bars 0.2–2 ms
tool.selectPCLasso([x1 y1; x2 y2; ...]);    % optional refine in PC space
tool.selectByLine([-0.1 250], [0.1 250]);   % select all waveforms crossing a line
newClusterId = tool.splitSelected();
tool.mergeWith(otherClusterId);
app.setThresholds([-80 250]);               % amplitude-window highlight, µV
```

## PETHs (main app → Tools → PETH / event alignment)

Open the **PETH** tool to see peri-event time histograms of the selected
cluster(s) aligned to behavioural events. Events load automatically from an
`events.mat` in the Phy folder, of the form:

```matlab
events.saccade = [...];   % alignment times in RECORDING seconds (same clock as spike_times/fs)
events.target  = [...];
events.reward  = [...];
events_cond.saccade = ["2","6",...];   % optional per-occurrence condition labels
meta = struct(...);                    % optional
```

The widget shows a grid of histogram PETHs — one row per selected cluster, one
column per chosen event (pick events in the list; set bin size and pre/post
window in ms). **Split by condition** replaces the bars with a smoothed spike
density function (Gaussian, mean ± SEM shading) per condition label (e.g. saccade
direction), which reads better than overlaid bars. Change the cluster selection
in the main window and click **Refresh**. A long `events.tsv`/`.csv` (columns
`event`, `time_s`, optional `condition`) is also accepted via **Load events…**.

## Dataset browser (main app → Tools → Dataset browser)

Open the **Dataset browser** to work through a whole dataset instead of one folder
at a time. Launch it two ways: from an open app (**Tools → Dataset browser**), or
**standalone** by running `DatasetBrowser` — the first time you open an entry
it launches a SpikeVisualizationApp for it and reuses that window afterwards. It
reads a dataset table (default `D:\CB paper\0-data\dataset.tsv`, one row per unit)
and shows it two ways — **Units** (one row per unit) or
**Recordings** (aggregated per `recording_tag`, with a unit count). Filter by
region / task / **Classification** (`bombcell_label`) / **Cell type**
(`c4_celltype`), type in the **Search** box (matches `recording_tag` and `unit`),
and click a column header to sort. Use **Columns…** to choose which columns the
Units view shows (e.g. add `in_original`); the choice is remembered per dataset.

**Open selected** (or double-click a row) loads that entry into the main app:

- a **recording** loads `<sortings root>\<recording_tag>\phy_postmerge`;
- a **unit** loads the sorting, selects that unit's cluster (resolved from its
  `merged_label` via `final_labels.tsv`) and opens the CurationTool on it — the
  main cluster table is still there to add other units of the recording.

On first use the browser asks you to confirm the **sortings root** (defaulting to
`…\curation_all`) and remembers it; change it later with **Sortings root…**. The
row currently loaded in the app is marked with a blue **▶**; rows with no sorting
folder are greyed out. A tag that has no Phy folder but does have an online (REX)
session falls back to loading that session spikes-only (see below).

**Curation labels & new units.** Click **Refresh** to re-read the saved sortings:
the `curation_label` and `curation_note` columns (sortable) are filled from each
recording's `cluster_curation.csv` — the SU/MU/good labels and notes you set in
SpikeVisualizationApp — and any cluster you labelled that isn't a dataset unit is
added as a new row (`unit = c<id>`, carrying only its recording_tag/region/task
plus the label/note; QC columns stay blank). The `cluster_id` column shows each
unit's cluster id, matching the app's Clusters table. (Refresh is manual so
opening the browser stays fast; your `dataset.tsv` is never changed.)

**Review status.** Each entry (unit or recording) carries a review state shown as
the row colour and in the `review` column: **unverified** (grey) · **WIP**
(yellow) · **issue** (red) · **discard** (dark) · **done** (green). Set it with the
**Set review status** buttons (and add a free-text **note**) for the selected row;
opening an entry auto-marks it WIP, and closing an app you opened from the browser
stamps it modified and recolours the row. State is saved to a **sidecar file**
next to the dataset (`<dataset>.review.tsv`) and merged back on load — your
original `dataset.tsv` is never modified, so it is always the revert point. Tick
**Lock** to make the browser read-only (no state is written).

## Spikes-only (online-sorted) sources

Some units were sorted **online** (REX) and have **only spike times** — no raw
trace and no waveforms. Open such a session with **Load file (.mat)…** and pick a
`<tag>_REX.mat` (or let the Dataset browser do it). The app detects the missing
waveforms/trace and hides the waveform, mean, PC-feature and raw-trace panels,
reflowing to show the **ISI** and **autocorrelogram** (and amplitude drift when
available) — everything that works from spike times alone. The same happens for
any Phy folder whose `data.dat` is absent.

## Saving

**Save sorting + labels** writes all per-spike arrays (`spike_times`,
`spike_clusters`, `spike_templates`, `amplitudes`, `pc_features`) back in one
time-sorted order so they stay aligned even after a realign, plus the labels for
every cluster to `cluster_curation.csv`. Writing is atomic (temp files then move,
with rollback) and retries through transient file locks; originals are backed up
once as `*.bak.npy`. On the next load the labels are read back from
`cluster_curation.csv`, so a reopened window matches exactly what you saved.
