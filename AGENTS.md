# AGENTS.md

## Purpose

This repository is a standalone R-side exporter for building KaroSpace-like HTML viewers from `.rds` inputs.

It is intentionally separate from the Python `KaroSpace` repository. The goal is to support R-native objects and workflows without modifying or depending on the `KaroSpace` repo at runtime.

## Hard Boundaries

- Do not modify the `KaroSpace` repository as part of work in this repo.
- Do not assume this repo can import Python code from `KaroSpace`.
- The viewer shell in this repo is a vendored snapshot, not a live dependency.
- Keep this repository independently usable and independently testable.

## Current Scope

Current implemented scope:

- Export from R-side objects to standalone HTML
- Supported inputs:
  - plain list inputs with `obs` and `coordinates`
  - `Seurat`
  - `SingleCellExperiment`
  - `SpatialExperiment`
- Render parity first:
  - section layout
  - categorical coloring
  - UMAP when present
  - core viewer interactions

Not fully implemented yet:

- robust gene-expression parity
- marker genes
- neighbor graph analytics
- interaction/contact markers
- desktop GUI

## Repo Structure

- `R/`
  Core exporter logic.
- `scripts/example_export.R`
  Fastest way to test a real `.rds` file.
- `scripts/karospace_build_r.R`
  CLI-style export entrypoint.
- `scripts/sync_viewer_shell.py`
  Refreshes the vendored viewer shell from the reference KaroSpace repo.
- `inst/viewer/karospace_viewer_shell.html`
  Vendored viewer shell used for HTML export.
- `tests/smoke_export.R`
  Minimal smoke test.

## Working Rules For Agents

- Prefer improving the existing R export path over introducing parallel code paths.
- Keep HTML export logic in `R/export.R` and payload construction in `R/payload.R`.
- Keep input normalization logic in `R/source.R`.
- If adding support for a new object type, normalize it into the same internal structure used by existing inputs.
- Preserve standalone export behavior: output should remain a self-contained HTML file.
- If syncing the viewer shell, use `scripts/sync_viewer_shell.py` rather than manually editing the vendored HTML unless the change is repo-local on purpose.
- Do not add a desktop app until the CLI and payload contract are stable.

## Seurat Notes

- Seurat spatial objects may contain multiple image slots in `@images`.
- Coordinate extraction must support multi-image integrated objects.
- Do not assume the first image covers all cells.
- Prefer normalized data when selecting an expression source, but avoid silent assumptions that hide missing data.

## Testing

Before finishing changes, run what is relevant:

```bash
Rscript tests/smoke_export.R
```

For real-data testing:

```bash
Rscript scripts/example_export.R --input /path/to/object.rds --inspect
Rscript scripts/example_export.R --input /path/to/object.rds
```

If auto-detection is wrong, rerun with explicit flags such as:

```bash
Rscript scripts/example_export.R \
  --input /path/to/object.rds \
  --groupby orig.ident \
  --initial-color seurat_clusters \
  --additional-colors integrated_snn_res.0.47
```

## Near-Term Priorities

Priority order for future work:

1. Gene-expression support and assay selection
2. Better metadata/default detection
3. Neighbor graph support
4. Marker-gene and interaction analytics
5. Performance improvements for large exports
6. Desktop application layer

## Output Expectations

When making changes, optimize for:

- correct export from real `.rds` objects
- stable payload structure
- viewer behavior matching the current vendored shell
- minimal surprise for users running the example script from Terminal
