## Pipeline stages

```
0. Build per-sample channels   - single sample, pre-assembled, or samplesheet
1. MYLOASM                     - long-read assembly (skipped if --assembly given)
2. MAPPING                     - minimap2 + samtools sort/index
3. BINNING                     - metabat2 / semibin2 (configurable models) / remag
4. DEREPLICATION                - dRep (metabat + semibin only, NOT remag)
5. CHECKM2_QC                   - bacterial completeness/contamination (NOT remag)
6. BUSCO_QC                     - eukaryotic QC (ALL bin sets, including remag)
7. GTDBTK_CLASSIFY               - bacterial taxonomy (NOT remag)          [toggle]
8. GVCLASS_TAXONOMY               - eukaryotic taxonomy (ALL bin sets)     [toggle]
   fast mode (SSU_EXTRACT_ALL or CMSEARCH_EUK) - alternative to 7 & 8      [toggle]
```

Steps 7, 8, and the fast-mode step are independently toggleable — run any
combination of them.

## Requirements

- Nextflow ≥ 23.10
- conda or mamba on every machine that will execute tasks — no
  Apptainer/Singularity/Docker required. Nextflow builds and caches a
  dedicated conda environment per process automatically, based on each
  process's own `conda` directive.
- `git` available on the host for the two pixi-managed tools (GVClass,
  SSUextract), which are cloned and built once, then cached via `storeDir`.

## Repo layout

```
main.nf                 workflow wiring, in the same order as the original .def file
nextflow.config          params, conda settings, resource labels
modules/
  myloasm.nf              assembly (params.read_type: hifi/ont preset)
  mapping.nf               minimap2 + samtools sort/index (params.read_type preset)
  binning.nf               DEPTH, METABAT2, SEMIBIN2, REMAG + BINNING subworkflow
  drep.nf                  CHECKM1_DB, DREP + DEREPLICATION subworkflow
  checkm2.nf               CHECKM2_DB, CHECKM2 + CHECKM2_QC subworkflow
  busco.nf                 UNZIP_BINS, BUSCO + BUSCO_QC subworkflow
  gtdbtk.nf                GTDBTK_DB, GTDBTK + GTDBTK_CLASSIFY subworkflow
  gvclass.nf               GVCLASS_SETUP, FNA_PREP (bin subsampling), COLLECT_SUBSAMPLE_REPORT, GVCLASS + GVCLASS_TAXONOMY subworkflow
  cmsearch.nf              fast-mode taxonomy: SSUEXTRACT_SETUP/SSU_EXTRACT (default) + CM_MODEL/CMSEARCH (legacy) subworkflows
```

## Key parameters

| Param | Default | Purpose |
|---|---|---|
| `input_fgz` | `/home/user/data/input.fastq.gz` | reads for single-sample mode |
| `nametag` | `""` | override the derived sample name (single-sample mode only) |
| `assembly` | `""` | path to an existing assembly — skips `MYLOASM` (single-sample mode only) |
| `samplesheet` | `""` | CSV path for multi-sample mode (columns: `fastq`, `assembly`) — see below. Takes priority over `input_fgz`/`assembly`/`nametag` when set |
| `read_type` | `'hifi'` | `'hifi'` or `'ont'` — selects the matching `MYLOASM`/`MAPPING` presets. See [Read type](#read-type-pacbio-hifi-vs-ont) below |
| `outdir` | `/home/user/data/output` | where results are published |
| `db_dir` | `'refdbs'` | **persistent, shared** location for one-time reference database/tool downloads (CheckM1, CheckM2, GTDB-Tk, GVClass, SSUextract) — deliberately separate from `outdir` so every run/sample reuses the same data instead of re-downloading. Relative to wherever you launch `nextflow` from — see the callout below |
| `run_metabat2` / `run_semibin2` / `run_remag` | `true` | toggle each binner independently |
| `semibin_models` | `'soil,self,global'` | comma-separated subset of SemiBin2 models to run |
| `run_gtdbtk` / `run_gvclass` | `true` | toggle full taxonomic classification steps |
| `run_cmsearch` | `false` | toggle the fast-mode alternative to GTDB-Tk/GVClass |
| `fast_mode_tool` | `'ssu_extract'` | `'ssu_extract'` (default, full taxonomy via SSUextract) or `'cmsearch'` (legacy, eukaryote-only, exact `fast_mode.def` parity). See [Fast mode](#fast-mode-ssuextract-vs-legacy-cmsearch) below |
| `ssuextract_db_profile` | `'curated'` | `'curated'` (~345 MiB, SILVA+PR2) or `'img'` (~841 MiB, + IMG sequences) |
| `ssuextract_tree_classification` | `false` | enable SSUextract's optional tree-based taxonomy refinement |
| `gvclass_max_bin_size_mb` | `50` | bins larger than this are **subsampled** (split into smaller parts), not excluded — see [Bin subsampling](#gvclass-bin-subsampling) below |

Sample names are always derived via `.simpleName` on the fastq filename
(strips all extensions) unless overridden by `--nametag` in
single-sample mode. Samplesheet mode has no per-row name override —
it's always `.simpleName` on that row's `fastq` path.

## Reference database location (`db_dir`)

`db_dir` defaults to `'refdbs'` — a relative path matching the `refdbs/`
folder already committed in this repo (kept empty via `.gitignore`; the
actual downloaded databases never get committed, only the folder
itself). This is a deliberate convenience default, with one thing worth
knowing: **a relative path resolves relative to wherever you run
`nextflow run` from, not relative to the repo.** As long as you always
launch the pipeline from the repo root (e.g. `cd` into it first, as
every example below assumes), `refdbs` reliably means the same folder
every time, and every run/sample shares the same cached databases.

If you ever launch this pipeline from a different working directory
(a script, a cron job, a different clone), pass an absolute path instead
to guarantee the same location regardless of where you start from:
```bash
--db_dir /home/you/genome_pipeline/Genome-Pipeline/refdbs
```

## Samplesheet format (multi-sample mode)

```csv
fastq,assembly
/data/sample01.fastq.gz,
/data/sample02.fastq.gz,/data/sample02_hifiasm_assembly.fa
/data/sample03.fastq.gz,
```

Leave `assembly` blank for any row that needs `MYLOASM` run; fill it in
for rows that already have one. All other pipeline settings (binner
toggles, thread counts, `db_dir`, `read_type`, fast-mode settings, etc.)
apply globally across every sample in the sheet.

## Read type: PacBio HiFi vs. ONT

`params.read_type` (`'hifi'` or `'ont'`, default `'hifi'`) controls the only
two places in the pipeline that are technology-specific:

| Step | HiFi | ONT |
|---|---|---|
| `MYLOASM` preset | `--hifi` | `--nano-r10` |
| `MAPPING` (minimap2 `-x`) | `map-hifi` | `map-ont` |

Everything downstream of `MAPPING` works purely off the assembly and/or
the resulting BAM coverage, with no awareness of which sequencing
technology produced them. An invalid value fails the run immediately
with a clear error rather than silently defaulting.

```bash
nextflow run main.nf --input_fgz sample.fastq.gz --outdir out/ --read_type ont -resume
```

**Note:** this is a single, pipeline-wide setting, including in
samplesheet mode. A mixed-technology batch (some HiFi, some ONT, in one
run) isn't currently supported — would need a `read_type` column added
to the samplesheet, similar to how `assembly` is already optional per row.

## Fast mode: SSUextract vs. legacy cmsearch

When `run_cmsearch` is enabled, `fast_mode_tool` selects which
implementation runs as the (much faster) alternative to full GTDB-Tk +
GVClass taxonomy:

- **`ssu_extract`** (default) — runs
  [SSUextract](https://github.com/NeLLi-team/ssuextract), a full
  Nextflow pipeline in its own right (same pixi-managed install pattern
  as GVClass). It detects 16S rRNA genes (routing through an RF00177-based
  index) and 18S rRNA genes (RF01960-based index) via bundled Infernal
  covariance models, then assigns taxonomy via BLAST against SILVA + PR2
  (`curated` profile) or SILVA + PR2 + IMG (`img` profile), with an
  optional tree-based classification refinement
  (`ssuextract_tree_classification`). Archaea are covered through SILVA's
  16S taxonomy rather than a separate archaeal model. Output lands in
  `<outdir>/<root>_ssuextract/results/`, including `cmsearch_summary.tsv`
  (the main per-hit taxonomy table) and `extracted/*.fna` (extracted SSU
  sequences).
- **`cmsearch`** (legacy) — the original `fast_mode.def` behavior exactly:
  a single cmsearch pass against the eukaryotic SSU rRNA model (Rfam
  RF01960) only, with hits extracted via `esl-sfetch`/`samtools faidx`
  into `<outdir>/<root>_extracted_ssu.fa`. No taxonomy assignment, no
  bacteria/archaea coverage — kept around for exact parity with the
  original script.

```bash
# fast mode, new default (SSUextract, full taxonomy, all 3 domains)
--run_cmsearch true

# fast mode, legacy single-domain cmsearch (exact fast_mode.def parity)
--run_cmsearch true --fast_mode_tool cmsearch
```

## GVClass bin subsampling

Bins larger than `gvclass_max_bin_size_mb` (GVClass has a hard input-size
limit and fails outright on oversized input) are **subsampled, not
excluded**: `FNA_PREP` splits them into multiple smaller files using
`gt splitfasta -targetsize` (from GenomeTools), and every resulting part
is still sent through to GVClass. A pipeline-wide
`<outdir>/gvclass_subsampled_bins_summary.tsv` records which bins were
split, their original size, and how many parts they became.

One edge case can't be fully resolved this way: a single contig larger
than the size limit can't be split further without cutting a sequence in
half, which `gt splitfasta` never does. That part is still passed through
to GVClass (rather than silently dropped), but flagged with a `[WARNING:
...]` note in the summary.

## Usage examples

**1. Basic single-sample run — full pipeline, default settings**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01 \
  --db_dir    refdbs \
  -resume
```

**2. Single sample with an explicit name**
```bash
nextflow run main.nf \
  --input_fgz /data/raw_reads_batch3_042.fastq.gz \
  --nametag   patient_A \
  --outdir    /results/patient_A \
  --db_dir    refdbs \
  -resume
```

**3. Pre-assembled input — skip MYLOASM**
```bash
nextflow run main.nf \
  --input_fgz /data/sample02.fastq.gz \
  --assembly  /data/sample02_hifiasm_assembly.fa \
  --outdir    /results/sample02 \
  --db_dir    refdbs \
  -resume
```
Still needs the reads (for mapping/depth calc), but `MYLOASM` never runs.

**4. Oxford Nanopore reads instead of PacBio HiFi**
```bash
nextflow run main.nf \
  --input_fgz /data/ont_sample.fastq.gz \
  --read_type ont \
  --outdir    /results/ont_sample \
  --db_dir    refdbs \
  -resume
```

**5. Fast mode — skip GTDB-Tk/GVClass, run SSUextract instead**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01_fast \
  --db_dir    refdbs \
  --run_gtdbtk   false \
  --run_gvclass  false \
  --run_cmsearch true \
  -resume
```

**6. Fast mode, legacy single-domain cmsearch**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01_fast_legacy \
  --db_dir    refdbs \
  --run_gtdbtk     false \
  --run_gvclass    false \
  --run_cmsearch   true \
  --fast_mode_tool cmsearch \
  -resume
```

**7. Only a subset of binners, and only some SemiBin2 models**
```bash
nextflow run main.nf \
  --input_fgz /data/sample01.fastq.gz \
  --outdir    /results/sample01 \
  --db_dir    refdbs \
  --run_remag       false \
  --semibin_models  soil,global \
  -resume
```

**8. Multi-sample batch run via samplesheet**
```bash
nextflow run main.nf \
  --samplesheet /data/batch1_samples.csv \
  --outdir      /results/batch1 \
  --db_dir      refdbs \
  -resume
```

**9. Combining samplesheet mode with fast mode and a binner subset**
```bash
nextflow run main.nf \
  --samplesheet     /data/batch1_samples.csv \
  --outdir          /results/batch1_fast \
  --db_dir          refdbs \
  --run_gtdbtk      false \
  --run_gvclass     false \
  --run_cmsearch    true \
  --semibin_models  global \
  -resume
```

**10. Re-running after a partial failure**
```bash
# some tasks failed partway through - just re-run the exact same command;
# -resume skips everything that already completed successfully and only
# retries what failed or never ran
nextflow run main.nf --samplesheet /data/batch1_samples.csv --outdir /results/batch1 --db_dir refdbs -resume
```

## Notes

- **`-resume` is a single dash**, not `--resume` — it's a Nextflow
  engine flag, not a pipeline param, so it doesn't go through `params.*`.
- **`db_dir` must point somewhere you have write access to** — a bare
  top-level path like `/refdbs` requires root to create and will fail.
  See [Reference database location](#reference-database-location-db_dir)
  above for the relative-vs-absolute tradeoff.
- **`db_dir` should point at the same location across every run** you
  want to share cached databases between — that's what makes the
  one-time downloads (CheckM1, CheckM2, GTDB-Tk, GVClass, SSUextract)
  actually pay off across samples/projects instead of re-downloading.
  These use Nextflow's `storeDir` directive, which checks whether the
  file already exists on disk (not just whether the task hash matches),
  so it skips re-downloading even on a brand-new run with no prior
  `work/` directory. The `refdbs/` folder in this repo is `.gitignore`'d
  so the (potentially many-GB) downloaded contents never get committed —
  only the empty folder itself is tracked.
- **dRep's internal dependency on checkM1** (a separate, older tool from
  CheckM2) needs its own reference database, downloaded once via
  `CHECKM1_DB` and registered per-task via `checkm data setRoot`. Without
  this, dRep fails at its filtering step with `!!! checkM failed !!!`
  and silently produces zero dereplicated genomes despite exiting status 0.
- Every binner toggle, taxonomy toggle, `semibin_models`, `read_type`,
  and fast-mode settings work identically in single-sample or
  samplesheet mode, since they're read once from `params.*` and applied
  uniformly regardless of how many samples are flowing through the channels.
- `GVCLASS_SETUP` and `SSUEXTRACT_SETUP` both clone a git repo and build
  a pixi-managed environment (not a conda package) into `db_dir` once,
  the same `storeDir`-cached way as the other databases.