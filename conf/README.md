# `conf/` — Execution and HPC configuration

This directory holds supplemental Nextflow configuration, separate from the primary `nextflow.config`. The intended `hpc.config` file defines optional executor profiles such as SLURM, PBS Pro, SGE, and LSF. Other configurations that are better handled in their own files can also be added here later.

- **Edit here:** To add or edit scheduler profiles; Add other configuration files related to nextflow orchestration. Currently we keep shared tool CPU and memory requirements in `nextflow.config`, but eventually they could be moved here.

For institutional examples see https://nf-co.re/configs and https://github.com/nf-core/configs. A profile file only takes effect when included by the main configuration. Installation tasks use the `install` label's local executor override.
