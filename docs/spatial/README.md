# Spatial image processing from segmentation exports

Functions live in `R/`; workflow `07_compute_spatial_data.R`
coordinates image processing, with `scripts/run_spatial_from_raw_alpine.sh` supplying
the preserved container. Shared input/intermediate paths are centralized in
`cfg$paths$spatial` in `config/config.R`, using the same path convention as other workflows.
Measurement settings remain in
`config/spatial_measurement_config.R`. This separate config preserves the existing
`config/spatial_config.R` used by other analyses; its updated goal-table labels differ
from that older configuration. Scientific kernels and thresholds are unchanged.

The launcher resolves paths from its repository location, not your current
working directory or a hard-coded macOS home directory:

- Tissue/nerve segmentation exports: `data/raw/tumor_stroma_nerve_masks/`
- Vessel segmentation exports: `data/raw/vessel_masks/`
- Generated single-channel TIFFs: `data/processed/spatial_masks/<run_name>/`
- Cleaned TIFFs, previews, geometry, provenance and comparisons: `data/processed/spatial_runs/<run_name>/`
- Only the three verified output CSVs: `data/processed/spatial/`

## Obtain and place the spatial TIFF inputs

The imaging-data source cited in the main README is
[BioImage Archive S-BIAD4129](https://www.ebi.ac.uk/bioimage-archive/galleries/S-BIAD4129/).
Download and extract the **tissue/nerve segmentation-mask TIFFs** and the
**vessel (CD31) segmentation-mask TIFFs** from the study's segmentation exports.
The archive's exact downloadable filenames/folder hierarchy have not been
verified here; identify the two sets by their content and the TIFF contract below.
The `02_segmentation_masks/tissue_and_nerve_masks/` and
`02_segmentation_masks/vessel_masks/` names in the project data-hierarchy notes
are organizational guidance, not verified paths within that public accession.

Copy the TIFF files themselves (not the enclosing downloaded folder) into these
locations relative to your Git repository root:

| Downloaded segmentation exports | Copy TIFFs into |
|---|---|
| Tumor/stroma/glass classes plus nerve instance labels | `data/raw/tumor_stroma_nerve_masks/` |
| CD31 vessel instance labels | `data/raw/vessel_masks/` |

For the local checkout, these resolve to
`~/Documents/GitHub/NMIBC_innervation_/data/raw/tumor_stroma_nerve_masks/` and
`~/Documents/GitHub/NMIBC_innervation_/data/raw/vessel_masks/`.
Use the same relative layout wherever the repository is checked out.

From the repository root, after extracting the downloads, replace the two source
paths below with their actual locations:

```bash
mkdir -p data/raw/tumor_stroma_nerve_masks data/raw/vessel_masks
cp /path/to/extracted/tissue_and_nerve_masks/*_binary_seg_maps.tif \
  data/raw/tumor_stroma_nerve_masks/
cp /path/to/extracted/vessel_masks/*_binary_seg_maps.tif \
  data/raw/vessel_masks/
```

For each requested image, both destination folders must directly contain
`<image_id>_binary_seg_maps.tif`; preserve image IDs and filenames. Do not add
another nested folder layer. Tissue/nerve exports must retain both pages:
page 1 is the semantic tissue map and page 2 is the nerve-instance map. Vessel
instances are read from page 2 of the CD31 export. Do not substitute `.im3`
images, spectrally unmixed component images, color overlays, or previously
processed single-channel masks. Although stored in `data/raw/`, these inputs
are segmentation exports rather than raw microscopy acquisitions.

The launcher generates the single-channel spatial inputs automatically under
`data/processed/spatial_masks/<run_name>/`. The separate filtered phenotype CSV
is also required to map nerve instance IDs to subtypes; pass its path explicitly
as described below. The TIFFs alone do not encode those subtype assignments.

Each original TIFF is named `<image_id>_binary_seg_maps.tif`. Tissue/nerve TIFF
page 1 contains semantic labels (0 tumor, 1 stroma, 2 glass); page 2 contains
nerve instances. Vessel TIFF page 2 contains vessel instances. These are
segmentation exports, not unsegmented microscopy. Nerve subtype assignments
require the filtered phenotype CSV with `image_id`, `objID`, `phenotype`;
phenotype values are subcluster_1, _2, _4, _5, _6, _7. objID is the TIFF instance
ID. Use the verified `all_nerve_component_phenotypes.csv` from July 6, or an
identical copy; neither clustering nor nearest-centroid reassignment is performed.

## Run spatial image processing on Alpine

Upload these files directly into the existing repository layout, not inside a
PUBLICATION_CODE folder. Supply the two segmentation directories above, the
phenotype CSV and the three original updated-format goal CSVs in a separate
reference directory. Preserve the successful RStudio container, overlay and
external library; see [SYSTEM_REQUIREMENTS.md](SYSTEM_REQUIREMENTS.md).

```bash
acompile --ntasks=1
cd /projects/YOUR_USERNAME/NMIBC_innervation_
bash scripts/run_spatial_from_raw_alpine.sh \
  /path/to/updated_format_goal_spatial_output \
  /path/to/all_nerve_component_phenotypes.csv \
  data/metadata/spatial/goal173.txt \
  spatial_processing173
```

Replace the example paths with the uploaded Alpine paths. Skip acompile if
already allocated. Default container paths can be overridden with
SPATIAL_CONTAINER, SPATIAL_OVERLAY and SPATIAL_R_LIBRARY pointing to the
preserved successful environment. The script checks the runtime, prepares fresh
TIFFs, computes all three tables, compares them inside that same container, and
copies them to data/processed/spatial only after every comparison passes. No goal
measurements are read during mask generation or calculation.

```bash
cat data/processed/spatial_runs/spatial_processing173/validation/REFERENCE_MATCH.txt
cat data/processed/spatial_runs/spatial_processing173/validation/PUBLISHED.txt
```

Both should report success. Inspect table_summary.csv for matching row keys and
column order, and column_comparison.csv for zero differing_cells. Areas/counts
and missingness are exact; other numeric values use 1e-8 + 1e-10 * abs(reference).
COMPLETE.txt alone is not verification. All three tables are required, including
object distances and vessel morphology that were not covered by the earlier
45-image band probe. The user reported an end-to-end PASS on five images. The organizational routing
change has local orchestration checks but has not itself been rerun on Alpine.

The default image manifest contains all 173 reference images. To process a
subset, supply a text file with one image ID per line and choose a new run name.
A subset run publishes only those images. Preserve
existing three spatial CSVs outside data/processed/spatial before publishing another
run; the launcher refuses overwrite. Other tables are unaffected.

Download the three verified CSVs to your local repository's data/processed/spatial/
and retain the corresponding run directory as provenance. Check downloaded
hashes against validation/published_tables_md5.csv (macOS: md5 -q FILE).

The included small data/metadata/spatial files record image selection and
supplied manual/clinical annotations, not measured outcomes. Original TIFFs,
phenotype data and reference tables are not included in this code update.
The separate existing scripts/08_run_spatial_analysis.R still reads its archived
RDS; this update does not alter its models or redirect its inputs.

The active launcher is `scripts/run_spatial_from_raw_alpine.sh`. Standalone
support scripts, small-test image lists and development provenance are archived
under `IGNORE/spatial_support/`; they are not required by the public workflow.

These are processed measurement inputs and run-provenance artifacts, not final manuscript statistics. Workflow 08 writes final model tables and figures to the main config’s spatial result root. Preprocessing-run previews remain alongside masks/geometry for QC; numbered workflow 07 runs the complete pipeline.

## Alpine runtime correction and validation status (2026-09-21)

The preserved Alpine container reports mmand **1.6.3**. The earlier hard check
for 1.7.0 came from the manuscript/local environment and was incorrect for this
container. Do not upgrade mmand to satisfy the old check. Any earlier description
of this version as locally tested is superseded by this observed Alpine version.
PASS was reported for five images from original segmentation exports through
all three goal tables, including distances and morphology, using mmand 1.6.3.
Returned files have not yet been independently reviewed; full-cohort reproduction
has not been established by this five-image test.

The Alpine launcher creates a temporary directory within the repository and
bind it at /tmp and /projects/$USER/.rstudioserver/rstudio-4.4.1/tmp_data inside
the container. The second mapping is necessary because R's system Renviron
hard-codes that path. The overlay and installed packages remain unchanged.
Run from /scratch/alpine/$USER (not the space-limited projects directory), with
the code and inputs copied there. Manual TMPDIR/APPTAINER_BINDPATH workarounds
are no longer needed in a fresh shell. The automated launcher change has local
shell syntax checks; the equivalent manual bind fix was used in the successful
Alpine run. Preserve completed outputs outside scratch for long-term storage.

## Result versus processed-data locations

The three spatial measurement CSVs are inputs to downstream analyses and remain
in `cfg$paths$spatial$measurements` (`data/processed/spatial/`). Generated TIFFs,
geometry and checkpoints remain under the configured processed mask/run roots.
The small PNG masks alongside TIFFs are technical previews tied to those data,
not manuscript figures. Workflow 08 statistical tables and plots use
`cfg$paths$results$spatial` with `tables/` and `figures/`, respectively. Historical border-overlay output locations are retained in the archived entry point.
No generated datasets, raw inputs or pre-existing output files were moved.

Two scientific configurations remain intentionally distinct: legacy
`config/spatial_config.R` and the five-image image processing's
`config/spatial_measurement_config.R`. Their display labels differ (including
`200-300 um` versus `201-300 um` and nerve subtype labels). They have not been
silently consolidated. Workflow 07 now runs the tested pipeline exclusively. The prior dual-route
entry point is preserved under `IGNORE/spatial_support/`. Workflow 08's
archived-input contract remains unchanged.

## Precomputed spatial measurements

`cfg$use_precomputed_spatial` in `config/config.R` defaults to `TRUE`. Workflow 07 then
skips the image-processing pipeline before loading its functions or image
packages, regardless of supplied command-line arguments.
Workflow 08 continues to read the precomputed dataset configured in
`config/spatial_config.R`; this switch does not run workflow 08 automatically.
Set `cfg$use_precomputed_spatial <- FALSE` explicitly before using the Alpine
image-processing launcher or running a fresh TIFF-to-table validation. Otherwise
workflow 07 reports that it skipped processing and generates no new measurements.

## Default environment for workflow 07

Run from the repository root:

```bash
Rscript scripts/07_compute_spatial_data.R GOALS PHENOTYPES [IMAGE_IDS] [RUN_NAME]
```

Replace the placeholders with paths (omit optional bracketed arguments when not
needed). With `cfg$use_precomputed_spatial = TRUE`, workflow 07 skips processing
without launching a container. With `FALSE`, it automatically invokes the Alpine
launcher; the container child performs the image processing and runtime checks.
The launcher sets an internal marker to prevent recursive container launches.
Do not set that marker manually. There is no fallback to local image processing.

A host R installation is needed to start this R entry point, but image-processing
packages are loaded inside the container. An Alpine compute allocation and the
preserved container, overlay and library must already be available. This does
not install packages, create an allocation or connect remotely. Direct use of
the shell launcher remains supported when host R is unavailable.
