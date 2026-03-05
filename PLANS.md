# KaroSpaceR – Implementation Plan

Goal: Achieve feature parity between the R export pipeline and the Python pipeline while ensuring the viewer supports gene expression, metadata defaults, spatial graphs, and scalable payload sizes.

---

## Phase 1: Gene Expression Parity

### Tasks
- [ ] Add a proper assay/layer resolver for each supported input type.
- [ ] For Seurat, default to normalized expression in `SCT` when available.
- [ ] Implement fallback order: `SCT` → `Spatial` → `counts`.
- [ ] Populate `available_genes` in the payload.
- [ ] Populate `genes_meta` in the payload.
- [ ] Export per-section gene matrices.
- [ ] Add CLI/example-script option `--assay`.
- [ ] Add CLI/example-script option `--genes`.

### Acceptance Criteria
- `Heart_A90_karospace_buildr.html` can switch between cluster coloring and true gene-expression coloring.

---

## Phase 2: Better Metadata and Viewer Defaults

### Tasks
- [ ] Improve auto-detection of `groupby`.
- [ ] Improve auto-detection of `initial_color`.
- [ ] Improve auto-detection of `additional_colors`.
- [ ] Improve detection of section-level metadata.
- [ ] Add optional `outline_by` support.
- [ ] Add cleaner metadata filtering for single-section datasets.
- [ ] Add cleaner metadata filtering for multi-section datasets.
- [ ] Normalize Seurat inputs to the standard payload structure.
- [ ] Normalize SingleCellExperiment inputs to the standard payload structure.
- [ ] Normalize SpatialExperiment inputs to the standard payload structure.
- [ ] Normalize plain-list inputs to the standard payload structure.

### Acceptance Criteria
- The same object exports with sensible defaults and minimal manual flags.

---

## Phase 3: Spatial Graph and Neighbor Features

### Tasks
- [ ] Add optional neighbor graph export.
- [ ] For Seurat, detect existing graph objects.
- [ ] Prefer existing graphs when available.
- [ ] Add option to derive graphs from spatial coordinates when missing.
- [ ] Populate `has_neighbors` in the payload.
- [ ] Export per-section edge lists.
- [ ] Export `neighbor_stats`.

### Acceptance Criteria
- Neighbor overlays work in the viewer.
- Neighbor summary panels function in the R export.

---

## Phase 4: Marker Genes and Interaction Features

### Tasks
- [ ] Implement marker-gene computation for categorical groupings.
- [ ] Add marker export to viewer payload.
- [ ] Add optional interaction/contact-style marker computation.
- [ ] Ensure interaction markers depend on the graph layer.
- [ ] Keep marker computations optional due to computational cost.

### Acceptance Criteria
- Viewer marker panels show real results instead of empty states.

---

## Phase 5: Performance and Payload Size

### Tasks
- [ ] Implement sparse gene encoding for expression matrices.
- [ ] Add optional array packing for large payloads.
- [ ] Add base64 encoding for large section payloads when needed.
- [ ] Implement a `lightweight export` mode.
- [ ] Reduce HTML payload size where possible.

### Acceptance Criteria
- Large datasets export successfully.
- HTML files remain usable and do not become excessively large.

---

## Status

- Phase 1: ⬜ Not started  
- Phase 2: ⬜ Not started  
- Phase 3: ⬜ Not started  
- Phase 4: ⬜ Not started  
- Phase 5: ⬜ Not started
