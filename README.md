# Genome-Pipeline

A Nextflow DSL2 long-read metagenomics pipeline: assembly, single-binner
binning, dereplication, and QC/taxonomy, finishing in one compiled
per-bin spreadsheet.

## Pipeline stages

```
0. Build per-sample channels   - single sample, pre-assembled, or samplesheet
1. MYLOASM                     - long-read assembly (skipped if --assembly given)
2. MAPPING                     - minimap2 + samtools sort/index
3. BINNING                     - exactly ONE binner runs: metabat2, semibin2, or remag
   (remag only)                  REMAG's rejected (bacterial-looking) contigs get
                                  re-binned by a second, user-chosen binner; both
                                  sets of bins are merged before dereplication
4. DEREPLICATION                - dRep, on whatever BINNING produced
                                  [--drep_ignore_genome_quality skips CheckM entirely]
5. CHECKM2_QC                   - completeness/contamination/assembly stats
6. BUSCO_QC                     - eukaryotic completeness QC
7. GTDBTK_CLASSIFY               - bacterial/archaeal taxonomy             [toggle]
8. GVCLASS_TAXONOMY               - eukaryotic/viral taxonomy             [toggle]
   fast mode (SSU_EXTRACT_ALL or CMSEARCH_EUK) - alternative to 7 & 8     [toggle]
9. REPORT                      - compiles CheckM2 + GTDB-Tk + GVClass + SSUextract
                                  into one per-bin spreadsheet (<root>_bin_report.xlsx)
```

Steps 7, 8, and the fast-mode step are independently toggleable. Only one
binner runs per invocation — see [Binner selection](#binner-selection)
below. The pipeline prints its total run time when it finishes — see
[Runtime announcement](#runtime-announcement) below.

## Requirements

- Nextflow (26.04.6 or later recommended)
- conda/mamba on `PATH` — every process builds its own environment automatically
- enough disk under `db_dir` for reference databases (CheckM1, CheckM2, GTDB-Tk, GVClass, SSUextract — tens of GB combined)

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
  cmsearch.nf                SSUextract (default) or legacy cmsearch, fast-mode alternative to 7 & 8
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
| `run_gtdbtk` / `run_gvclass` | `true` | toggle full taxonomic classification steps |
| `run_cmsearch` | `false` | enable fast mode (alternative to GTDB-Tk/GVClass) |
| `fast_mode_tool` | `'ssu_extract'` | `'ssu_extract'` (default) or `'cmsearch'` (legacy) |
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

## Fast mode: SSUextract vs. legacy cmsearch

`--run_cmsearch true` enables a fast alternative to full taxonomic
classification (`GTDBTK_CLASSIFY`/`GVCLASS_TAXONOMY`), extracting and
classifying 16S/18S rRNA genes directly from the assembly instead.
`fast_mode_tool` picks which implementation:

- `'ssu_extract'` (default): the actively maintained
  [SSUextract](https://github.com/NeLLi-team/ssuextract) pipeline,
  covering bacteria, archaea, *and* eukaryotes (16S + 18S).
- `'cmsearch'` (legacy): a single Rfam covariance model (RF01960),
  eukaryotes only.

Fast mode operates directly on the assembly — it never depends on
`BINNING` or `DEREPLICATION`, so it runs independently of whichever
binner you chose.

```bash
--run_cmsearch true                          # default tool: ssu_extract
--run_cmsearch true --fast_mode_tool cmsearch  # legacy single-model path
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
- **SSUextract's summary** (if `run_cmsearch`) — `ssuextract_tax` and
  `SSU_seq`, attached to whichever bin actually contains the contig the
  rRNA hit came from (left blank for every other bin)

Any toggle that's off, or a tool output that's missing or has drifted in
format, fills that section's columns with `N/A` rather than failing the
whole report.

**Two things worth knowing about this report that aren't fully
battle-tested yet:**
- GVClass's exact taxonomy column name isn't hardcoded — the script
  searches for any column containing "taxonomy" or "classification" in
  GVClass's own header. Worth double-checking it picked the right one
  against a real run.
- The SSUextract-to-bin matching (scanning every bin fasta's headers for
  the rRNA hit's contig ID) is this pipeline's own heuristic, not
  something SSUextract or any other tool does natively — it reproduces
  the pattern of "one bin gets the hit, others are blank," but hasn't
  been validated against a real multi-bin run yet.

## Usage examples

**All paths below are relative**, matching this repo's own convention
(`data/` for inputs, `results/` for outputs, `refdbs/` for the shared
database cache). A relative path resolves relative to wherever you run
`nextflow run` from — these examples assume you're launching from the
repo root every time. Use an absolute path instead for anything that
needs to resolve the same way regardless of your current directory.

**1. Basic single-sample run — metabat2 (default), full taxonomy**
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

**5. Fast mode only, skip full taxonomy**
```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir    results/sample01_fast \
  --db_dir    refdbs \
  --run_gtdbtk     false \
  --run_gvclass    false \
  --run_cmsearch   true \
  -resume
```

**6. Legacy cmsearch fast mode**
```bash
nextflow run main.nf \
  --input_fgz data/sample01.fastq.gz \
  --outdir     results/sample01_fast_legacy \
  --db_dir     refdbs \
  --run_gtdbtk    false \
  --run_gvclass   false \
  --run_cmsearch  true \
  --fast_mode_tool cmsearch \
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

**8. Multi-sample, fast mode only**
```bash
nextflow run main.nf \
  --samplesheet      data/batch1_samples.csv \
  --outdir           results/batch1_fast \
  --db_dir           refdbs \
  --run_gtdbtk       false \
  --run_gvclass      false \
  --run_cmsearch     true \
  -resume
```

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
deeper problem, especially for long-running steps like `SSU_EXTRACT` or
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
