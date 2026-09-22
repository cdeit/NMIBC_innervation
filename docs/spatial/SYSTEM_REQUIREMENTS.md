# Spatial image-processing system requirements

## Validated environment

These are observed values from the successful CU Alpine RStudio container, not
suggestions to install the latest versions.

| Component | Recorded requirement |
|---|---|
| R | 4.4.1 |
| R platform | x86_64-pc-linux-gnu |
| Operating system inside container | Ubuntu 22.04.5 LTS |
| EBImage | 4.48.0 |
| fftwtools | 0.9-11 |
| FFTW runtime | 3.3.8; Ubuntu libfftw3-double3:amd64 3.3.8-2ubuntu8 |
| fftwtools linkage | Dynamic libfftw3.so.3 from the Ubuntu environment |
| tiff | 0.1-12 |
| png (mask previews) | 0.1-8 in the recorded Alpine session |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| tibble | 3.2.1 |
| CPU in the initial successful Alpine probe | AMD EPYC 7543 |
| Shell/container tools | Bash, Apptainer, an allocated Alpine compute session |
| Vessel morphology and optional preprocessing | mmand 1.6.3 (observed in the successful Alpine run) |

Object-border distances use EBImage and do not require RANN. The updated vessel
CSV also contains morphology, so mmand is required by this runner. The returned
band-only Alpine probe did not validate the morphology calculations. A five-image end-to-end PASS was reported with mmand 1.6.3, including
the vessel morphology and object-distance columns.
No Python, Arrow, glmmTMB, emmeans or plotting packages are needed.

## Existing Alpine container and library

The launcher uses these defaults from the successful probe:

```text
$CURC_CONTAINER_DIR_OOD/rstudio-server-4.4.1.sif
/projects/$USER/.rstudioserver/rstudio-4.4.1/rstudio-server-4.4.1_overlay.img
/projects/$USER/Rstudio_libs/4.4.1
```

It mounts the personal overlay read-only and sets R_LIBS_USER. Override
`SPATIAL_CONTAINER`, `SPATIAL_OVERLAY`, `SPATIAL_R_LIBRARY` only to locate the same
environment elsewhere. Required input/output paths are explicitly bound into the
container. Output parent directories must exist. The launcher performs no package
installation, updates, or remote login. An RStudio web session is not required.
Do not substitute a batch R module or Alpine Linux (an unrelated distribution).

The preprocessing and measurement scripts call `check_runtime()` in `R/runtime.R`
to verify core versions, R/platform, Ubuntu release, FFTW package and dynamic linkage. It fails
on a different runtime. Development-only direct execution can explicitly set
`SPATIAL_ALLOW_UNVALIDATED=1`; the Alpine launcher disables that override.

## Why package names alone are insufficient

FFT-based convolution can produce residuals near zero. The original cleanup
uses a strict computed `glass_contact > 0` rule. A numerical difference can change
whether an entire tissue component is retained and therefore alter border bands
and distances. The code intentionally preserves that rule. Installing the same R
package versions on another OS/CPU does not establish numerical agreement.

The observed CPU is provenance, not proof that this particular CPU is the sole
requirement. The exact cause was not isolated to one package, CPU feature or
library build. Preserve the successful environment and validate outputs.

## What to preserve with a publication

- Exact container image or an accessible immutable image identifier plus SHA256.
- Personal overlay AND external R package library snapshot: the container alone
  does not necessarily contain all loaded packages.
- SHA256 hashes of those preserved artifacts (do not invent a container digest).
- Per-run sessionInfo, package versions/binary hashes, ldd linkage, CPU and OS
  records, input hashes, source snapshot and resolved spatial configuration.
- Explicit image list and access to the source segmentation/phenotype data.

The runner automatically records the per-run items under `provenance/`. An renv
lockfile alone cannot preserve FFTW, the OS, or CPU/runtime behavior. Alternative container environments have not been validated.

Memory depends on image dimensions and the number of connected border segments:
the historical distance function retains one distance raster per segment. Use a
serial job and an allocation with enough memory; no universal memory/time minimum
has been validated. Do not infer whole-cohort runtime from the small band probe.

Official Alpine RStudio documentation:
https://curc.readthedocs.io/en/latest/open_ondemand/rstudio.html

## Output storage

The source-compatible signed TIFF writer is uncompressed. Preprocessing writes
12 full-resolution single-channel TIFFs per image; geometry export writes seven
more. Allow sufficient project/scratch storage for your image sizes and cohort.
PNG previews and RDS files are compressed. No GPU or graphical display is needed.

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
