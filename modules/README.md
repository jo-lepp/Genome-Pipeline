# `modules/` — Nextflow processes and workflows

This directory contains the DSL2 modules that implement the pipeline's analysis and installation steps. `main.nf` imports and these processes into the full workflow.

- **Edit here:** To change commands, staged inputs/outputs, tool-specific options, or analysis behavior.
- **Handled elsewhere:** CPU/memory defaults and executor selection belong in `nextflow.config` and `conf/hpc.config`.
- **Examples:** `install.nf` defines software/database installers; `gvclass.nf`, `cmsearch.nf`, and `gtdbtk.nf` define analysis steps.

When moving an installed tool or database, check both its installer in `install.nf` and the corresponding analysis module. Keep input/output declarations consistent with `main.nf`.
