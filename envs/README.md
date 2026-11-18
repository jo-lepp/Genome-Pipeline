# `envs/` — Conda environment recipes

Each YAML file specifies a Conda environment used by one or more Nextflow processes. Modules select the correct environment with their `conda` directive; Nextflow creates or reuses the corresponding environment in `conda_cache_dir`.

- **Edit here:** To add, remove, or pin packages needed by a process. Update package versions, or pin them down to a particular version.
- **Special case:** `pixi_tools.yml` supplies Git and Pixi; GVClass and SSUextract install their own nested environments under `software/`.

If you add a new environment, reference it from its module. To make it part of the explicit `--install_software` workflow, also add its name to the `environmentNames` list in `modules/install.nf`.
