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
| `SpikeVisualizationApp.m` | Main app: cluster list, waveforms, mean±CI, amplitude drift, ISI, ACG, PC features, and a scrollable raw-trace view with spike markers. Cluster labelling and save. |
| `ISISplitTool.m` | Companion window: lasso ISI bars near 0 → highlight those spikes in the PC-feature scatter → split them into a new cluster. |
| `LoadSpikesPhy.m` | Reader for a Phy export folder (`params.py`, `spike_*.npy`, `pc_features.npy`, `cluster_*.tsv`), extracting one waveform per spike from the raw data channel. |

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

## Splitting bad spikes out of a cluster (ISI → PC features)

1. Select the cluster in the list and click **ISI split...**.
2. **Circle ISI bars** – draw a freehand region over the bars close to 0 in the
   ISI histogram. Every spike taking part in an ISI inside that range is
   highlighted (orange) in the PC-feature scatter (temporal PCA of the ch0
   waveform — the same features Phy's FeatureView shows).
3. **Lasso PC features** – optionally draw a lasso around the offending branch in
   the scatter (red) to refine exactly which spikes to remove.
4. **Split off selected** – the selection is moved to a new cluster in the main
   app, where you can judge it (waveform, ISI, ACG) and keep or delete it.

The same actions are scriptable for reproducible curation:

```matlab
app  = SpikeVisualizationApp(phyDir);
tool = ISISplitTool(app, clusterId);
tool.selectIsiRange([0 2]);                 % ISI bars 0–2 ms
tool.selectPCLasso([x1 y1; x2 y2; ...]);    % optional refine in PC space
newClusterId = tool.splitSelected();
```

## Saving

**Save sorting + labels** writes the edited cluster assignment back to
`spike_clusters.npy` (backing up the original to `spike_clusters.bak.npy`) and
the SU/MU/noise labels to `cluster_curation.csv` in the Phy folder.
