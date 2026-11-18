# Genome-Pipeline

A Nextflow DSL2 long-read metagenomics pipeline: assembly, single-binner
binning, dereplication, and QC/taxonomy, finishing in one compiled
per-bin spreadsheet.

## Pipeline stages

```
0. Build per-sample channels    - single sample, pre-assembled, or samplesheet
1. MYLOASM                      - long-read assembly (skipped if --assembly given)
2. MAPPING                      - minimap2 + samtools sort/index
3. BINNING                      - exactly ONE binner runs: metabat2, semibin2, or remag
   (remag only)                  REMAG's rejected (bacterial-looking) contigs get
                                  re-binned by a second, user-chosen binner; both
                                  sets of bins are merged before dereplication
4. DEREPLICATION                - dRep, on whatever BINNING produced
                                  [--drep_ignore_genome_quality skips CheckM entirely]
5. CHECKM2_QC                   - completeness/contamination/assembly stats
6. BUSCO_QC                     - eukaryotic completeness QC
7. GTDBTK_CLASSIFY              - bacterial/archaeal taxonomy             [toggle]
8. GVCLASS_TAXONOMY             - eukaryotic/viral taxonomy               [toggle]
9. SSU screening (SSU_EXTRACT_ALL or CMSEARCH_EUK) - 16S/18S rRNA, runs   [toggle]
   alongside 7 & 8 by default - see Fast mode below
10. REPORT                     - compiles CheckM2 + GTDB-Tk + GVClass + SSU screening
                                  into one per-bin spreadsheet (<root>_bin_report.xlsx)
```

Steps 7, 8, and 9 are independently toggleable, and all three run by
default. Only one binner runs per invocation — see
[Binner selection](#binner-selection) below. The pipeline prints its
total run time when it finishes — see
[Runtime announcement](#runtime-announcement) below.

## Requirements

- Nextflow (26.04.6 or later recommended)
- conda/mamba on `PATH` — every process builds its own environment automatically
- enough disk under `db_dir` for reference databases (CheckM1, CheckM2, GTDB-Tk, GVClass, SSUextract — tens of GB combined)

## Installation

Genome-Pipeline provides automated installation of its software dependencies and reference databases through Nextflow. Software installation and database installation are managed separately so that users can reuse existing installations, share reference databases across projects, and avoid unnecessary downloads.

The pipeline uses **Conda for most software dependencies** and **Pixi for GVClass and SSUextract**. Singularity, Apptainer, and Docker are not required.

### Prerequisites

Before installing Genome-Pipeline, ensure that the following are available:

- Nextflow (version 26.04.6 or later recommended)
- Conda or Mamba, accessible from the command line
- Git and standard Unix utilities
- Internet access to download software packages and reference databases
- Sufficient disk space for reference databases, particularly GTDB-Tk and CheckM2

Clone the repository and enter the project directory:

```bash
git clone https://github.com/hhollandmoritz/Genome-Pipeline.git
cd Genome-Pipeline
```

If using the development installation branch:

```bash
git checkout installation
```

### Installation directories

Genome-Pipeline stores Conda environments, Pixi applications, and reference databases in separate locations. These can be configured in `nextflow.config`:

```groovy
params {
    conda_cache_dir = "${projectDir}/conda_cache"
    software_dir    = "${projectDir}/software"
    db_dir          = "${projectDir}/refdbs"
}
```

| Parameter | Purpose | Default location |
|---|---|---|
| `conda_cache_dir` | Conda environments managed by Nextflow | `conda_cache/` |
| `software_dir` | GVClass and SSUextract Pixi installations | `software/` |
| `db_dir` | Reference databases used by the pipeline | `refdbs/` |

These directories are persistent and can be reused between pipeline runs.

Users working on HPC systems may prefer absolute paths to shared storage locations. For example:

```bash
nextflow run main.nf \
  --install_all true \
  --conda_cache_dir /shared/software/genome_pipeline/conda_cache \
  --software_dir /shared/software/genome_pipeline/software \
  --db_dir /shared/databases/genome_pipeline
```

**Important:** Pixi environments are not relocatable. GVClass and SSUextract should remain in their original installation directories after installation. Moving these directories can cause Pixi to fail.

### Installing software and databases

The pipeline supports three installation modes.

| Parameter | Action |
|---|---|
| `--install_all true` | Install all software environments and reference databases |
| `--install_software true` | Install Conda environments and Pixi applications only |
| `--install_databases true` | Install reference databases using previously installed software |

**Option 1: Complete installation (recommended for new users)**

Install all dependencies and databases in one command:

```bash
nextflow run main.nf \
  --install_all true
```

Nextflow initializes the Conda environments, installs GVClass and SSUextract through Pixi, and downloads the reference databases once the relevant software is available.

**Option 2: Install software only**

```bash
nextflow run main.nf \
  --install_software true
```

This installs the pipeline's Conda environments and the Pixi-managed GVClass and SSUextract applications, without downloading their reference databases. It is useful if you want to create symlinks to a shared set of databases.

**Option 3: Install databases only**

```bash
nextflow run main.nf \
  --install_databases true
```

This installs the reference databases without repeating the software installation workflow.

The required Conda and Pixi installations must already be available. For a fresh installation, run `--install_software true` first or use `--install_all true`.


### Reusing existing installations

Users who already have software or reference databases installed can configure the pipeline to use their existing locations.

For example:

```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir results/sample01 \
  --db_dir /shared/databases/genome_pipeline \
  -resume
```
The `db_dir` parameter specifies the parent directory containing all reference databases.

Genome-Pipeline expects the following directory structure, which is created by `--install_databases true`:

```text
refdbs/
├── checkm1_db/
│   └── checkm_data/
│       └── ...                      # CheckM1 reference data
├── checkm2_data/
│   └── checkm2_db/
│       └── CheckM2_database/
│           └── uniref100.KO.1.dmnd
├── gtdbtk_db/
│   └── gtdbtk_data/
│       ├── metadata/
│       ├── markers/
│       ├── masks/
│       ├── msa/
│       ├── pplacer/
│       ├── radii/
│       └── taxonomy/
├── busco_data/
│   └── busco_datasets/
│       └── ...                      # Selected BUSCO lineage datasets
├── gvclass/
│   └── ...                          # GVClass reference database
├── ssuextract/
│   └── ...                          # SSUextract reference database
└── rfam_cm/
    └── RF01960.cm
```
If the system you are working on already has some or all of these databases installed and prepared for these tools somewhere else, you can use symlinks to link to them on the system. 

```
cd Genome-Pipeline
ln -s /path/to/gtdb-tk/dbs /path/to/db_dir/gtdbtk_db/
```
For shared installations, users should ensure that the required files are readable by the accounts running the analysis.


## Repo layout

```
main.nf                    entry workflow
nextflow.config             params, resource labels, conda settings
modules/
  myloasm.nf                MYLOASM
  mapping.nf                MAPPING
  binning.nf                DEPTH, METABAT2, SEMIBIN2, REMAG + BINNING subworkflow
                             (single-binner selection + REMAG reject rebinning)
  drep.nf                   CHECKM1_DB, DREP + DEREPLICATION subworkflow
                             (ignoreGenomeQuality support, bin-format normalization)
  checkm2.nf                CHECKM2_DB, CHECKM2 + CHECKM2_QC subworkflow
  busco.nf                  BUSCO + BUSCO_QC subworkflow
  gtdbtk.nf                 GTDBTK_DB, GTDBTK + GTDBTK_CLASSIFY subworkflow
  gvclass.nf                GVCLASS_SETUP, FNA_PREP, GVCLASS, COLLECT_SUBSAMPLE_REPORT + GVCLASS_TAXONOMY subworkflow
  cmsearch.nf                SSUextract (default) or legacy cmsearch, SSU screening alongside full taxonomy
  report.nf                  REPORT: compiles the final per-bin spreadsheet
```

## Key parameters

| Param | Default | Notes |
|---|---|---|
| `input_fgz` | `/home/user/data/input.fastq.gz` | reads for single-sample mode |
| `nametag` | `""` | sample name override (single-sample mode) |
| `assembly` | `""` | pre-existing assembly — skips MYLOASM if given |
| `samplesheet` | `""` | enables multi-sample mode — see below |
| `read_type` | `'hifi'` | `'hifi'` or `'ont'` — picks minimap2's preset |
| `outdir` | `/home/user/data/output` | where results are published |
| `db_dir` | `'refdbs'` | shared reference-database cache, relative to launch dir |
| `binner` | `'metabat'` | `'metabat'`, `'semibin'`, or `'remag'` — see [Binner selection](#binner-selection) |
| `remag_reject_binner` | `'metabat'` | `'metabat'` or `'semibin'` — only used when `binner` is `'remag'` |
| `semibin_models` | `'soil'` | comma-separated SemiBin2 model(s); used when `binner` (or `remag_reject_binner`) is `'semibin'` |
| `drep_ignore_genome_quality` | `false` | skip CheckM entirely in dRep — see [dRep genome quality filtering](#drep-genome-quality-filtering) |
| `run_gtdbtk` | `true` | bacterial/archaeal taxonomy |
| `run_gvclass` | `true` | eukaryotic/viral taxonomy |
| `run_ssu_screen` | `true` | 16S/18S rRNA screening — see [Fast mode](#fast-mode-ssu-screening) |
| `fast_mode_tool` | `'ssu_extract'` | `'ssu_extract'` (default) or `'cmsearch'` (legacy) — which tool `run_ssu_screen` uses |
| `fast_mode` | `false` | preset: switches to `cmsearch`, disables GVClass — see [Fast mode](#fast-mode-ssu-screening) |
| `gvclass_max_bin_size_mb` | `50` | bins larger than this are **subsampled**, not excluded — see [Bin subsampling](#gvclass-bin-subsampling) |
| `busco_lineage` | `'eukaryota_odb10'` | BUSCO lineage dataset |

## Reference database location (`db_dir`)

All one-time downloads (CheckM1, CheckM2, GTDB-Tk, GVClass, SSUextract)
are cached under `db_dir`, relative to wherever you launch `nextflow run`
from — see [Usage examples](#usage-examples) for the relative-path
convention this repo follows. CheckM1's database is only downloaded when
`drep_ignore_genome_quality` is `false` (the default) — if you set it to
`true`, that download is skipped entirely since dRep never calls CheckM.

## Samplesheet format (multi-sample mode)

```csv
fastq,assembly
data/sample01.fastq.gz,
data/sample02.fastq.gz,data/sample02_hifiasm_assembly.fa
data/sample03.fastq.gz,
```
Sample name = `fastq_file.simpleName`. No per-row name override in
samplesheet mode.

## Read type: PacBio HiFi vs. ONT

`--read_type hifi` (default) uses minimap2's `map-hifi` preset;
`--read_type ont` uses `map-ont`. This also affects MYLOASM's assembly
mode internally.

## Binner selection

Only **one** binner runs per invocation, controlled by `--binner`:

| `binner` value | What runs |
|---|---|
| `'metabat'` (default) | metabat2 only |
| `'semibin'` | SemiBin2, using whichever model(s) `semibin_models` lists |
| `'remag'` | REMAG (eukaryote-focused) — see [REMAG reject rebinning](#remag-reject-rebinning) below |

```bash
# default
--binner metabat

# run SemiBin2 instead, with two models
--binner semibin --semibin_models soil,global

# run REMAG, re-binning its bacterial rejects with SemiBin2
--binner remag --remag_reject_binner semibin
```

This is a deliberate simplification from earlier versions of this
pipeline, which ran all three binners in parallel and used DAS_Tool to
pick a consensus set. With one binner chosen up front, every downstream
step (`DEREPLICATION`, `CHECKM2_QC`, `BUSCO_QC`, `GTDBTK_CLASSIFY`,
`GVCLASS_TAXONOMY`) runs exactly once per sample, on exactly one bin set
— no more per-label branching.

## REMAG reject rebinning

REMAG's niche is eukaryotic bins — by design, it filters out
bacterial-looking contigs rather than binning them. When `binner` is
`'remag'`, those rejected contigs (written to `*_non_eukaryotic.fasta`
via REMAG's own `--save-filtered-contigs` flag) are **not discarded**:
they're re-binned by a second binner of your choice
(`remag_reject_binner`, `'metabat'` or `'semibin'`), reusing the same
BAM already computed for the sample — no re-alignment needed, since
`jgi_summarize_bam_contig_depths` only reports depth for whichever
contigs are present in the given fasta.

Both sets of bins — REMAG's own eukaryotic bins, and the secondary
binner's bacterial bins — are merged into one list before
`DEREPLICATION`. Because they represent different organisms, they won't
share meaningful ANI overlap, so dRep simply keeps both untouched rather
than trying to dereplicate them against each other.

```bash
--binner remag --remag_reject_binner metabat
```

## dRep genome quality filtering

By default, dRep runs CheckM internally and filters out genomes below
its completeness/contamination thresholds before clustering. CheckM's
models are trained on bacteria and archaea — they will often
**incorrectly score real eukaryotic genomes as low-quality and filter
them out**, which defeats the purpose of running REMAG (or any
eukaryote-focused workflow) in the first place.

Setting `--drep_ignore_genome_quality true` passes dRep's own
`--ignoreGenomeQuality` flag, which skips CheckM entirely — no
quality-based filtering happens, and genomes are clustered purely on
ANI. This also means the CheckM1 database is never downloaded for that
run (see [Reference database location](#reference-database-location-db_dir)
above).

```bash
# recommended whenever binner is 'remag', or whenever eukaryotic
# genomes are expected in the sample
--drep_ignore_genome_quality true
```

This is **off by default** — only enable it when you specifically expect
genomes CheckM would misjudge; for an all-bacterial sample, the default
quality filtering is still the safer choice.

## Fast mode: SSU screening

By default, this pipeline runs **everything**: full bacterial/archaeal
taxonomy (`GTDB-Tk`), full eukaryotic/viral taxonomy (`GVClass`), *and*
lightweight 16S/18S rRNA screening (`run_ssu_screen`) — the three aren't
mutually exclusive, and all three run side by side unless you turn one
off.

`run_ssu_screen` controls whether that rRNA screening step runs at all;
`fast_mode_tool` picks which implementation it uses:

- `'ssu_extract'` (default): the actively maintained
  [SSUextract](https://github.com/NeLLi-team/ssuextract) pipeline,
  covering bacteria, archaea, *and* eukaryotes (16S + 18S).
- `'cmsearch'` (legacy): a single Rfam covariance model (RF01960),
  eukaryotes only.

This screening step operates directly on the assembly — it never
depends on `BINNING` or `DEREPLICATION`, so it runs independently of
whichever binner you chose, and independently of GTDB-Tk/GVClass too.

```bash
--run_ssu_screen false   # turn screening off entirely
--fast_mode_tool cmsearch  # use the legacy single-model tool instead of SSUextract
```

### The `--fast_mode` preset

If you actually want a fast, lightweight run — screening only, no full
taxonomy — `--fast_mode true` is a single flag that does three things at
once, overriding whatever `fast_mode_tool`/`run_gvclass`/`run_ssu_screen`
were otherwise set to:

1. switches `fast_mode_tool` to `'cmsearch'` (the lighter tool)
2. disables `GVCLASS_TAXONOMY` entirely
3. forces `run_ssu_screen` on (so the fast path actually runs something)

`GTDB-Tk` is **not** affected by `fast_mode` — it keeps whatever
`run_gtdbtk` says, on by default. If you want a true minimal run, turn
it off explicitly too:

```bash
# fast preset: cmsearch screening, GVClass off, GTDB-Tk still on by default
--fast_mode true

# fast preset + no bacterial/archaeal taxonomy either
--fast_mode true --run_gtdbtk false
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

## Output spreadsheet

The final step, `REPORT`, compiles everything into one per-bin workbook:
`<outdir>/<root>_bin_report.xlsx`. It joins:

- **CheckM2's `quality_report.tsv`** — the backbone of the report:
  `Completeness`, `Contamination`, `Coding_Density`, `Contig_N50`,
  `Average_Gene_Length`, `Genome_Size`, `GC_Content`,
  `Total_Coding_Sequences`, `Total_Contigs`, `Max_Contig_Length`
- **GTDB-Tk's `summary.tsv`** (if `run_gtdbtk`) — `GTDB_classification`,
  plus `closest_genome_reference`/`_ani`/`_af` (GTDB-Tk's `fastani_*`
  columns, falling back to `closest_placement_*` when ANI-based
  classification wasn't available) and `red_value`
- **GVClass's `gvclass_summary.tsv`** (if `run_gvclass`) —
  `GVclass_classification`
- **SSU screening's summary** (if `run_ssu_screen`) — `ssuextract_tax`
  and `SSU_seq`, attached to whichever bin actually contains the contig
  the rRNA hit came from (left blank for every other bin)

Any toggle that's off, or a tool output that's missing or has drifted in
format, fills that section's columns with `N/A` rather than failing the
whole report.

**Two things worth knowing about this report that aren't fully
battle-tested yet:**
- GVClass's exact taxonomy column name isn't hardcoded — the script
  searches for any column containing "taxonomy" or "classification" in
  GVClass's own header. Worth double-checking it picked the right one
  against a real run.
- The SSU-screening-to-bin matching (scanning every bin fasta's headers
  for the rRNA hit's contig ID) is this pipeline's own heuristic, not
  something SSUextract/cmsearch or any other tool does natively — it
  reproduces the pattern of "one bin gets the hit, others are blank,"
  but hasn't been validated against a real multi-bin run yet.

## Usage examples

**All paths below are relative**, matching this repo's own convention
(`data/` for inputs, `results/` for outputs, `refdbs/` for the shared
database cache). A relative path resolves relative to wherever you run
`nextflow run` from — these examples assume you're launching from the
repo root every time. Use an absolute path instead for anything that
needs to resolve the same way regardless of your current directory.

**1. Basic single-sample run — everything on by default**

metabat2 binning, full GTDB-Tk + GVClass taxonomy, and SSUextract
screening all run:
```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir    results/sample01 \
  --db_dir    refdbs \
  -resume
```

**2. REMAG, with bacterial rejects re-binned by SemiBin2, ignoring genome quality in dRep**
```bash
nextflow run main.nf \
  --input_fgz data/raw_reads_batch3_042.fastq.gz \
  --outdir    results/patient_A \
  --db_dir    refdbs \
  --binner               remag \
  --remag_reject_binner  semibin \
  --drep_ignore_genome_quality true \
  -resume
```

**3. Pre-supplied assembly, SemiBin2 with two models**
```bash
nextflow run main.nf \
  --input_fgz data/sample02.fastq.gz \
  --assembly  data/sample02_hifiasm_assembly.fa \
  --outdir    results/sample02 \
  --db_dir    refdbs \
  --binner          semibin \
  --semibin_models  soil,global \
  -resume
```

**4. ONT reads**
```bash
nextflow run main.nf \
  --input_fgz data/ont_sample.fastq.gz \
  --outdir    results/ont_sample \
  --db_dir    refdbs \
  --read_type ont \
  -resume
```

**5. Fast preset — screening only, no full taxonomy at all**
```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir    results/sample01_fast \
  --db_dir    refdbs \
  --fast_mode true \
  --run_gtdbtk false \
  -resume
```

**6. Full taxonomy only, screening turned off**
```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir     results/sample01_full_tax_only \
  --db_dir     refdbs \
  --run_ssu_screen false \
  -resume
```

**7. Multi-sample run via samplesheet**
```bash
nextflow run main.nf \
  --samplesheet     data/batch1_samples.csv \
  --outdir          results/batch1 \
  --db_dir          refdbs \
  --binner          metabat \
  -resume
```

**8. Multi-sample, fast preset**
```bash
nextflow run main.nf \
  --samplesheet      data/batch1_samples.csv \
  --outdir           results/batch1_fast \
  --db_dir           refdbs \
  --fast_mode        true \
  --run_gtdbtk       false \
  -resume
```
## Working on HPC systems

The pipeline supports execution on high-performance computing (HPC) systems through Nextflow's built-in job scheduler integration. Nextflow automatically submits individual pipeline processes as separate jobs, using the CPU and memory requirements specified in `nextflow.config`.

By default the pipeline runs in "local" mode(no job scheduler). Users running on HPC systems can select and configure their job scheduler using a Nextflow profile.

### Selecting a job scheduler

Scheduler profiles are defined in `conf/hpc.config` and loaded automatically by `nextflow.config`.

| Scheduler | Profile | Command option |
|---|---|---|
| Local (default) | None required | No `-profile` needed |
| SLURM | `slurm` | `-profile slurm` |
| PBS Pro | `pbs` | `-profile pbs` |
| Sun Grid Engine | `sge` | `-profile sge` |
| LSF | `lsf` | `-profile lsf` |

For example, to run the pipeline using SLURM:

```bash
nextflow run main.nf \
  -profile slurm \
  --input_fgz data/sample01.fastq.gz \
  --outdir results/sample01 \
  -resume
```

Installation tasks are local-only due to many worker-nodes being unable to connect to the internet. If your system allows more flexibility, you can relax this requirement in the by commenting out the lines below in the 'process' block of `nextflow.config`. 
```
withLabel: 'install' {
        executor = 'local'
    }
```

### Institution-specific HPC configurations

HPC systems often require additional settings, such as job partitions, allocation accounts, submission limits, or cluster-specific resource requirements.

The [nf-core institutional configuration repository](https://nf-co.re/configs) provides an extensive collection of example Nextflow configurations for HPC systems worldwide. The corresponding configuration files are available at [github.com/nf-core/configs](https://github.com/nf-core/configs).

Users can either manually edit the `conf/hpc.config` file or adapt an appropriate configuration for their institution and supply it using Nextflow's `-c` option:

```bash
nextflow run main.nf \
  -profile slurm \
  -c my_cluster.config \
  --input_fgz data/sample01.fastq.gz \
  --outdir results/sample01 \
  -resume
```

For example, a SLURM cluster configuration might contain:

```groovy
process {
    queue = 'normal'
    clusterOptions = '--account=my_account'
}

executor {
    queueSize = 20
}
```

Extra configuration files can be supplied:

```bash
nextflow run main.nf \
  -profile slurm \
  -c my_cluster.config,my_resources.config \
  --input_fgz data/sample01.fastq.gz \
  --outdir results/sample01 \
  -resume
```

The `-c` option merges user configurations with the pipeline defaults. Existing parameters and resource requirements from nextflow.config will be used unless overridden in the config files supplied by the user.

### Installing software and databases on HPC systems

The pipeline separates software and reference database installation from analysis.

Installation processes use the local executor, regardless of whether or not an HPC profile is specified. This allows users to perform installation on login or build nodes where software installation is permitted.

Once installation is complete, subsequent analysis runs do not invoke installation processes unless installation is explicitly requested.

Conda environments, Pixi-managed applications, and reference databases must be accessible from the compute nodes. Users can configure their locations through `conda_cache_dir`, `software_dir`, and `db_dir` in `nextflow.config` or through pipeline parameters.

### Launching Nextflow through the scheduler

On some HPC systems, it is preferable or required to launch the main Nextflow process through a job that manages the submission of all the other jobs rather than keeping it running on a login node.

For SLURM, users can create a batch script that starts Nextflow with `-profile slurm`. The main batch job manages the workflow, while Nextflow submits individual analysis processes as additional SLURM jobs. 

The main Nextflow job must remain running until the workflow completes, but can be launched with minimal compute requirements.

For additional information, consult the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html) and the [nf-core HPC configuration examples](https://nf-co.re/configs).


## Runtime announcement

Every run prints a short summary when it finishes, success or failure,
via Nextflow's `workflow.onComplete` handler:

```
Pipeline finished at : 2026-10-09T17:42:11.123-04:00
Total run time        : 4h 32m 10s
Success               : true
Exit status           : 0
```

`Total run time` is Nextflow's own `workflow.duration` — wall-clock time
for the whole invocation, not CPU time summed across tasks. On a
`-resume`'d run, this reflects only the time this particular invocation
actually spent running (cached tasks don't add to it), not the combined
time across every attempt that led up to a finished result.

Under Nextflow's strict parser, this handler must be declared as an
**assignment** inside the entry workflow — `workflow.onComplete = { ... }`
— not the older `workflow.onComplete { ... }` method-call form, which
compiles but misbehaves at runtime.

## Notes

### Troubleshooting: SSH/VPN drop mid-run

A dash `-` exit status with an empty `.command.err` usually means the
outer Nextflow process itself was killed — an SSH disconnect, VPN
drop, VM reboot, or manually killing the wrong process.

**Confirm this is what happened:**
```bash
ps aux | grep nextflow
```
If the only match is the `grep` command matching its own search term
(no real `nextflow run` process listed), the pipeline isn't running at
all anymore — that's the confirmation.

**The fix is simply to re-launch with `-resume`:**
```bash
nextflow run main.nf --input_fgz ... --outdir ... --db_dir refdbs -resume
```
Everything that completed successfully before the crash is reused from
cache. The one task that was actually running when the crash happened
will need to re-run from scratch (Nextflow's caching works at the
whole-task level, not mid-task) — this is expected, not a sign of a
deeper problem, especially for long-running steps like SSU screening or
`GVCLASS` where getting partway through before a crash is likely.

### Troubleshooting: GTDB-Tk / pplacer memory errors

If `GTDBTK` fails with a pplacer exception, this is a documented
GTDB-Tk issue: pplacer's own multi-threading can make the OS report far
more memory in use than is actually needed, which can trigger an
out-of-memory kill. GTDB-Tk's own FAQ recommends `--scratch_dir` and
`--pplacer_cpus 1` — both are already set in this pipeline's `GTDBTK`
process. If it still happens, check how many CPUs are available on the
host versus what `threads_gtdbtk` requests.

### Troubleshooting: AVHE VPN required

`as-amoliver-vm2.syr.edu` is unreachable without the Syracuse VPN. Use
`tmux` to persist sessions across disconnects:
```bash
tmux new -s pipeline       # start a session
tmux attach -t pipeline    # reattach after reconnecting
# inside tmux: Ctrl+b c (new window), Ctrl+b n/p (switch), Ctrl+b d (detach)
```

### Troubleshooting: `work/` directory size

Nextflow's `work/` directory accumulates every task's files, including
cached ones from old runs. Use `nextflow clean -f -before <run_name>` to
prune old runs safely without touching ones you might still want to
`-resume`.
