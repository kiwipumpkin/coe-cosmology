> ⚠️ DRAFT — work in progress.

# figure_sources

- `src/` — scripts and input data needed by the packet figures that were not already elsewhere in the repo, mirroring their original location (`/workspace/<path>` → `src/<path>`).
- `plotdata/figNNN_<image>.npz` — the data points actually plotted, captured when each figure was saved (one array per line/scatter/image/bar set).
- `PATHMAP.json` — original path → location in this repo (relative to `draft/current/`), for every script and input file.
- `OUTDIRS.json` — output directories the scripts write to.
- `rebuild_layout.py <root>` — recreates the original layout under `<root>` from the repo so any figure script can be rerun.

See [`../FIGURE_MANIFEST.md`](../FIGURE_MANIFEST.md).
