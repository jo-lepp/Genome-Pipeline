# Genome-Pipeline

information on pipeline here (come back to this)

- Note that the first iteration of this pipeline was build for apptainer. .def files are found in the apptainer_defs folder

## Pipeline stages

```
0. Build per-sample channels   - single sample, pre-assembled, or samplesheet
1. MYLOASM                     - long-read assembly (skipped if --assembly given)
2. MAPPING                     - minimap2 + samtools sort/index
3. BINNING                     - metabat2 / semibin2 (configurable models) / remag
4. DEREPLICATION                - dRep (metabat + semibin only, NOT remag)
5. CHECKM2_QC                   - bacterial completeness/contamination (NOT remag)
6. BUSCO_QC                     - eukaryotic QC (ALL bin sets, including remag)
7. GTDBTK_CLASSIFY              - bacterial taxonomy (NOT remag)          [toggle]
8. GVCLASS_TAXONOMY              - eukaryotic taxonomy (ALL bin sets)     [toggle]
   CMSEARCH_EUK                  - fast-mode alternative to 7 & 8         [toggle]
```

Steps 7, 8, and `CMSEARCH_EUK` are independently toggleable — run any
combination of them.

## Requirements

- Nextflow ≥ 23.10
- conda or mamba on every machine that will execute tasks — no
  Apptainer/Singularity/Docker required. Nextflow builds and caches a
  dedicated conda environment per process automatically, based on each
  process's own `conda` directive.

## Repo layout

```
main.nf                 workflow wiring, in the same order as the original .def file
nextflow.config          params, conda settings, resource labels
modules/
  myloasm.nf              assembly
  mapping.nf              minimap2 + samtools sort/index
  binning.nf               DEPTH, METABAT2, SEMIBIN2, REMAG + BINNING subworkflow
  drep.nf                  DREP + DEREPLICATION subworkflow
  checkm2.nf               CHECKM2_DB, CHECKM2 + CHECKM2_QC subworkflow
  busco.nf                 UNZIP_BINS, BUSCO + BUSCO_QC subworkflow
  gtdbtk.nf                GTDBTK_DB, GTDBTK + GTDBTK_CLASSIFY subworkflow
  gvclass.nf               GVCLASS_SETUP, FNA_PREP, COLLECT_SKIPPED, GVCLASS + GVCLASS_TAXONOMY subworkflow
  cmsearch.nf              CM_MODEL, CMSEARCH + CMSEARCH_EUK subworkflow (fast mode)
```

## Key parameters

| Param | Default | Purpose |
|---|---|---|
| `input_fgz` | `/home/user/data/input.fastq.gz` | reads for single-sample mode |
| `nametag` | `""` | override the derived sample name (single-sample mode only) |
| `assembly` | `""` | path to an existing assembly — skips `MYLOASM` (single-sample mode only) |
| `samplesheet` | `""` | CSV path for multi-sample mode — see below. Takes priority over `input_fgz`/`assembly`/`nametag` when set |
| `outdir` | `/home/user/data/output` | where results are published |
| `db_dir` | `/refdbs` | **persistent, shared** location for one-time reference database downloads (CheckM2, GTDB-Tk, GVClass) — deliberately separate from `outdir` so every run/sample reuses the same databases instead of re-downloading |
| `run_metabat2` / `run_semibin2` / `run_remag` | `true` | toggle each binner independently |
| `semibin_models` | `'soil,self,global'` | comma-separated subset of SemiBin2 models to run |
| `run_gtdbtk` / `run_gvclass` | `true` | toggle full taxonomic classification steps |
| `run_cmsearch` | `false` | toggle the fast-mode alternative (Rfam SSU rRNA search) |
| `gvclass_max_bin_size_mb` | `50` | bins larger than this are excluded from GVClass input (it fails outright on oversized bins) |
| `read_type` | `'hifi'` | `'hifi'` or `'ont'` — selects the matching `MYLOASM`/`MAPPING` presets. See [Read type: PacBio HiFi vs. ONT](#read-type-pacbio-hifi-vs-ont) below |

Sample names are always derived via `.simpleName` on the fastq filename
(strips all extensions) unless overridden by `--nametag` in
single-sample mode. Samplesheet mode has no per-row name override —
it's always `.simpleName` on that row's `fastq` path.

## Samplesheet format (multi-sample mode)

```csv
fastq,assembly
/data/sample01.fastq.gz,
/data/sample02.fastq.gz,/data/sample02_hifiasm_assembly.fa
/data/sample03.fastq.gz,
```

Leave `assembly` blank for any row that needs `MYLOASM` run; fill it in
for rows that already have one. All other pipeline settings (binner
toggles, thread counts, `db_dir`, etc.) apply globally across every
sample in the sheet.

## Read type: PacBio HiFi vs. ONT

`params.read_type` (`'hifi'` or `'ont'`, default `'hifi'`) controls the only
two places in the pipeline that are technology-specific:

| Step | HiFi | ONT |
|---|---|---|
| `MYLOASM` preset | `--hifi` | `--nano-r10` |
| `MAPPING` (minimap2 `-x`) | `map-hifi` | `map-ont` |

Everything downstream of `MAPPING` — `BINNING`, `DEREPLICATION`,
`CHECKM2_QC`, `BUSCO_QC`, `GTDBTK_CLASSIFY`, `GVCLASS_TAXONOMY`, and
`CMSEARCH_EUK` — works purely off the assembly and/or the resulting BAM
coverage, with no awareness of which sequencing technology produced
them, so none of those steps need to change based on this setting.
(`SEMIBIN2` in particular already runs with the generic
`--sequencing-type long_read` flag, not anything HiFi-specific.)

An invalid value (anything other than `'hifi'`/`'ont'`) fails the run
immediately with a clear error, rather than silently falling back to
one of the two presets.

```bash
# PacBio HiFi (default - no flag needed)
nextflow run main.nf --input_fgz sample.fastq.gz --outdir out/ --db_dir refdbs/

# Oxford Nanopore R10
nextflow run main.nf --input_fgz sample.fastq.gz --outdir out/ --db_dir refdbs/ --read_type ont
```

**Note:** this is a single, pipeline-wide setting — including in
samplesheet mode, where it applies to every sample in the batch. If you
need a *mixed* batch (some samples HiFi, some ONT, in one run), the
samplesheet format would need a `read_type` column added per row,
similar to how the `assembly` column is already optional per row. Not
currently supported, but a small extension if you need it.

## Usage examples

**1. Basic single-sample run — full pipeline, default settings**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01 \
  --db_dir    /refdbs \
  -resume
```
Runs `MYLOASM` → `MAPPING` → all three binners → dRep → CheckM2 → BUSCO
→ GTDB-Tk → GVClass. `root` is derived from the filename (`sample01`),
since `--nametag` isn't set.

**2. Single sample with an explicit name**
```bash
nextflow run main.nf \
  --input_fgz /data/raw_reads_batch3_042.fastq.gz \
  --nametag   patient_A \
  --outdir    /results/patient_A \
  --db_dir    /refdbs \
  -resume
```
Output folders are named `patient_A_myloasm/`, `patient_A_metabat/`,
etc., instead of the messy filename-derived name.

**3. Pre-assembled input — skip MYLOASM**
```bash
nextflow run main.nf \
  --input_fgz /data/sample02.fastq.gz \
  --assembly  /data/sample02_hifiasm_assembly.fa \
  --outdir    /results/sample02 \
  --db_dir    /refdbs \
  -resume
```
Still needs the reads (for mapping/depth calc), but `MYLOASM` never
runs — the assembly you supply is used directly downstream.

**4. Fast mode — skip GTDB-Tk/GVClass, run cmsearch instead**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01_fast \
  --db_dir    /refdbs \
  --run_gtdbtk   false \
  --run_gvclass  false \
  --run_cmsearch true \
  -resume
```

**5. Only a subset of binners, and only some SemiBin2 models**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01 \
  --db_dir    /refdbs \
  --run_remag       false \
  --semibin_models  soil,global \
  -resume
```
Runs metabat2 + SemiBin2 (soil and global models only, skipping
self-supervised) — REMAG never runs at all, so BUSCO/GVClass only see 3
bin sets instead of 5.

**6. Multi-sample batch run via samplesheet**
```bash
nextflow run main.nf \
  --samplesheet /data/batch1_samples.csv \
  --outdir      /results/batch1 \
  --db_dir      /refdbs \
  -resume
```

**7. Combining samplesheet mode with fast mode and a binner subset**
```bash
nextflow run main.nf \
  --samplesheet     /data/batch1_samples.csv \
  --outdir          /results/batch1_fast \
  --db_dir          /refdbs \
  --run_gtdbtk      false \
  --run_gvclass     false \
  --run_cmsearch    true \
  --semibin_models  global \
  -resume
```
All the toggles are global across the whole batch — every sample in
the CSV gets the same settings.

**8. Re-running after a partial failure**
```bash
# some tasks failed partway through (e.g. one sample OOM'd on CheckM2) - just
# re-run the exact same command; -resume skips everything that already
# completed successfully and only retries what failed or never ran
nextflow run main.nf --samplesheet /data/batch1_samples.csv --outdir /results/batch1 --db_dir /refdbs -resume
```

## Notes

- **`-resume` is a single dash**, not `--resume` — it's a Nextflow
  engine flag, not a pipeline param, so it doesn't go through `params.*`.
- **`--db_dir` should point at the same path across every run** you
  want to share cached databases (CheckM2, GTDB-Tk, GVClass) between —
  that's what makes the one-time downloads actually pay off across
  samples/projects instead of re-downloading. These use Nextflow's
  `storeDir` directive, which checks whether the file already exists on
  disk (not just whether the task hash matches), so it skips
  re-downloading even on a brand-new run with no prior `work/` directory.
- Every binner toggle (`run_metabat2`, `run_semibin2`, `run_remag`),
  taxonomy toggle (`run_gtdbtk`, `run_gvclass`, `run_cmsearch`), and
  `semibin_models` work identically in single-sample or samplesheet
  mode, since they're read once from `params.*` and applied uniformly
  regardless of how many samples are flowing through the channels.
- `gvclass_max_bin_size_mb` before handing it to GVClass, and writes a
  pipeline-wide `<outdir>/gvclass_skipped_bins_summary.tsv` listing
  exactly which bins were excluded and why.


