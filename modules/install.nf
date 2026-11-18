
/*
 * Genome-Pipeline installation workflows (Nextflow DSL2)
 *
 * INSTALL_ENVS      : Initialize environments from envs/*.yml
 * INSTALL_DATABASES : Install reference databases and Pixi projects
 * INSTALL_ALL       : Initialize environments, then install databases
 *
 * Conda environments are managed by Nextflow and stored in
 * conda.cacheDir, configured in nextflow.config.
 *
 * Reference data and Pixi project checkouts are stored under
 * params.db_dir.
 *
 * No Singularity or Apptainer containers are required.
 */

nextflow.enable.dsl = 2


// ============================================================================
// CONDA ENVIRONMENT INSTALLATION
// ============================================================================

/*
 * Nextflow creates the environment automatically before the process runs.
 * The environment is cached using conda.cacheDir.
 *
 * Each task accepts a YAML recipe and reports successful initialization.
 */

process INSTALL_ENV {
    tag "$name"
    label 'install'

    conda { env_yaml }

    input:
    tuple val(name), val(env_yaml)

    output:
    val(name), emit: ready

    script:
    """
    set -euo pipefail
    echo "Conda environment initialized: ${name}"
    echo "Recipe: ${env_yaml}"
    echo "Environment prefix: \$CONDA_PREFIX"
    """
}

// ============================================================================
// SOFTWARE INSTALLATION
// ============================================================================
// for GVClass and SSUextract, which are nested Pixi projects

process INSTALL_GVCLASS_SOFTWARE {
    label 'install'
    conda "${projectDir}/envs/pixi_tools.yml"
    storeDir "${params.software_dir}/.install_state"

    input:
    val ready

    output:
    path 'gvclass_setup_complete', emit: ready

    script:
    """
    set -euo pipefail

    INSTALL_DIR="${params.software_dir}/gvclass"
    mkdir -p "${params.software_dir}"

    if [ ! -d "\$INSTALL_DIR/.git" ]; then
        git clone '${params.gvclass_repo_url}' "\$INSTALL_DIR"
    fi

    cd "\$INSTALL_DIR"
    pixi install --frozen

    touch "\$OLDPWD/gvclass_setup_complete"
    """
}

process INSTALL_SSUEXTRACT_SOFTWARE {
    label 'install'
    conda "${projectDir}/envs/pixi_tools.yml"
    storeDir "${params.software_dir}/.install_state"

    input:
    val ready

    output:
    path 'ssuextract_setup_complete', emit: ready

    script:
    """
    set -euo pipefail

    INSTALL_DIR="${params.software_dir}/ssuextract"
    mkdir -p "${params.software_dir}"

    if [ ! -d "\$INSTALL_DIR/.git" ]; then
        git clone '${params.ssuextract_repo_url}' "\$INSTALL_DIR"
    fi

    cd "\$INSTALL_DIR"

    sed -i '/^[[:space:]]*time[[:space:]]*=/d' config/base.config

    pixi install --frozen

    touch "\$OLDPWD/ssuextract_setup_complete"
    """
}

// ============================================================================
// REFERENCE DATABASE INSTALLATION
// ============================================================================

// ----------------------------------------------------------------------------
// CheckM1 (called by drep - not default)
// ----------------------------------------------------------------------------

process INSTALL_CHECKM1_DB {
    label 'install'
    storeDir "${params.db_dir}/checkm1_db"

    conda "${projectDir}/envs/drep.yml"

    input:
    val ready

    output:
    path 'checkm_data', emit: db

    script:
    """
    set -euo pipefail

    mkdir -p checkm_data
    wget -O checkm_data.tar.gz '${params.checkm1_data_url}'

    tar -xzf checkm_data.tar.gz -C checkm_data
    rm checkm_data.tar.gz
    """
}


// ----------------------------------------------------------------------------
// CheckM2
// ----------------------------------------------------------------------------

process INSTALL_CHECKM2_DB {
    label 'install'
    storeDir "${params.db_dir}/checkm2_data"

    conda "${projectDir}/envs/checkm2.yml"

    input:
    val ready

    output:
    path 'checkm2_db', emit: db

    script:
    """
    set -euo pipefail

    mkdir -p checkm2_db

    checkm2 database --download \
        --path "\$(pwd)/checkm2_db"
    """
}


// ----------------------------------------------------------------------------
// GTDB-Tk
// ----------------------------------------------------------------------------

process INSTALL_GTDBTK_DB {
    label 'install'
    storeDir "${params.db_dir}/gtdbtk_db"

    conda "${projectDir}/envs/gtdbtk.yml"

    input:
    val ready

    output:
    path 'gtdbtk_data', emit: db

    script:
    """
    set -euo pipefail

    mkdir -p gtdbtk_data

    wget -O gtdbtk_data.tar.gz \
        '${params.gtdbtk_data_url}'

    tar -xzf gtdbtk_data.tar.gz \
        -C gtdbtk_data \
        --strip-components=1

    rm gtdbtk_data.tar.gz
    """
}


// ----------------------------------------------------------------------------
// BUSCO
// ----------------------------------------------------------------------------

process INSTALL_BUSCO_DB {
    label 'install'
    storeDir "${params.db_dir}/busco_data"

    conda "${projectDir}/envs/busco.yml"

    input:
    val ready

    output:
    path 'busco_datasets', emit: db

    script:
    """
    set -euo pipefail

    mkdir -p busco_datasets

    busco --download \
        '${params.busco_lineage ?: 'eukaryota_odb10'}' \
        --download_path "\$(pwd)/busco_datasets"
    """
}


// ----------------------------------------------------------------------------
// GVClass
// ----------------------------------------------------------------------------

/*
 * Conda provides Git and Pixi.
 *
 * The GVClass project is cloned into the Nextflow task directory.
 * Pixi installs its own locked environment inside the checkout.
 *
 * setup-db downloads reference resources into the project checkout,
 * which is persisted under params.db_dir.
 */

process INSTALL_GVCLASS_DB {
    label 'install'
    conda "${projectDir}/envs/pixi_tools.yml"
    storeDir "${params.db_dir}/gvclass"

    input:
    val ready

    output:
    path 'gvclass_db_complete', emit: db

    script:
    """
    set -euo pipefail

    INSTALL_DIR="${params.software_dir}/gvclass"
    DB_DIR="${params.db_dir}/gvclass"

    if [ ! -d "\$INSTALL_DIR/.pixi/envs/default" ]; then
        echo "ERROR: GVClass software is not installed." >&2
        echo "Run: nextflow run main.nf --install_software true" >&2
        exit 1
    fi

    mkdir -p "\$DB_DIR"
    cd "\$INSTALL_DIR"

    pixi run setup-db --database "\$DB_DIR"

    touch "\$OLDPWD/gvclass_db_complete"
    """
}


// ----------------------------------------------------------------------------
// SSUextract
// ----------------------------------------------------------------------------

/*
 * SSUextract is a nested Nextflow/Pixi project.
 *
 * Its repository and downloaded resources are stored in params.db_dir.
 *
 * The nested Nextflow task time limits are removed.
 */

process INSTALL_SSUEXTRACT_DB {
    label 'install'
    conda "${projectDir}/envs/pixi_tools.yml"
    storeDir "${params.db_dir}/ssuextract"

    input:
    val ready

    output:
    path 'ssuextract_db_complete', emit: db

    script:
    """
    set -euo pipefail

    INSTALL_DIR="${params.software_dir}/ssuextract"
    DB_DIR="${params.db_dir}/ssuextract"

    if [ ! -d "\$INSTALL_DIR/.pixi/envs/default" ]; then
        echo "ERROR: SSUextract software is not installed." >&2
        echo "Run: nextflow run main.nf --install_software true" >&2
        exit 1
    fi

    mkdir -p "\$DB_DIR"
    cd "\$INSTALL_DIR"

    pixi run setup \
        --database_profile '${params.ssuextract_db_profile ?: 'curated'}' \
        --database_path "\$DB_DIR"

    touch "\$OLDPWD/ssuextract_db_complete"
    """
}


// ----------------------------------------------------------------------------
// Rfam RF01960 covariance model
// ----------------------------------------------------------------------------

process INSTALL_RFAM_DB {
    label 'install'
    storeDir "${params.db_dir}/rfam_cm"

    input:
    val ready

    output:
    path 'RF01960.cm', emit: db

    script:
    """
    set -euo pipefail

    wget -O RF01960.cm \
        '${params.rfam_eukarya_cm_url}'

    test -s RF01960.cm
    """
}


// ============================================================================
// WORKFLOW: INSTALL_ENVS
// ============================================================================

workflow INSTALL_ENVS {

    main:

    /*
     * All Conda environment definitions used by the pipeline.
     * Additional YAML files can be added to this inventory.
     */

    def environmentNames = [
        'myloasm',
        'mapping',
        'metabat2',
        'semibin',
        'remag',
        'drep',
        'checkm2',
        'busco',
        'gtdbtk',
        'genometools',
        'infernal',
        'dastool',
        'pixi_tools',
        'report'
    ]

    // Verify that every YAML recipe exists.
    environmentNames.each { name ->
        def recipe = file("${projectDir}/envs/${name}.yml")
        if (!recipe.exists()) {
            error("Missing Conda environment definition: ${recipe}")
        }
    }

    // Supply absolute recipe paths to Nextflow.
    ch_envs = Channel.fromList(
        environmentNames.collect { name ->
            tuple(
                name,
                file("${projectDir}/envs/${name}.yml")
                    .toAbsolutePath()
                    .toString()
            )
        }
    )

    INSTALL_ENV(ch_envs)

    // install software pixi environments for gvclass and ssuextract
    // depends on pixi environment being installed
    ch_pixi_ready = INSTALL_ENV.out.ready
        .filter { it == 'pixi_tools' }
        .first()

    INSTALL_GVCLASS_SOFTWARE(ch_pixi_ready)
    INSTALL_SSUEXTRACT_SOFTWARE(ch_pixi_ready)

        emit:
        ready = INSTALL_ENV.out.ready
        software_ready = INSTALL_GVCLASS_SOFTWARE.out.ready
            .map { 'gvclass' }
            .mix(
                INSTALL_SSUEXTRACT_SOFTWARE.out.ready
                    .map { 'ssuextract' }
            )
    }


// ============================================================================
// WORKFLOW: INSTALL_DATABASES
// ============================================================================

workflow INSTALL_DATABASES {

    take:
    ch_ready

    main:

    // When a container is ready, it produces ch_ready.
    // Databases are installed when ch_ready is available.

    ch_drep = ch_ready
        .filter { it == 'drep' }
        .first()

    ch_checkm2 = ch_ready
        .filter { it == 'checkm2' }
        .first()

    ch_gtdbtk = ch_ready
        .filter { it == 'gtdbtk' }
        .first()

    ch_busco = ch_ready
        .filter { it == 'busco' }
        .first()

    ch_pixi = ch_ready
        .filter { it == 'pixi_tools' }
        .first()

    ch_gvclass = ch_ready
        .filter { it == 'gvclass' }
        .first()

    ch_ssuextract = ch_ready
        .filter { it == 'ssuextract' }
        .first()
    
    ch_infernal = ch_ready
        .filter { it == 'infernal' }
        .first()

    // Each database process starts independently as soon as
    // its corresponding environment is initialized.

    INSTALL_CHECKM1_DB(ch_drep)
    INSTALL_CHECKM2_DB(ch_checkm2)
    INSTALL_GTDBTK_DB(ch_gtdbtk)
    INSTALL_BUSCO_DB(ch_busco)

    // These depend on the Pixi conda environment AND the 
    // software installation in the software_dir to run.
    // Pixi projects don't like to have moving directory structures
    // so the software folder is the fix.
    INSTALL_GVCLASS_DB(ch_gvclass)
    INSTALL_SSUEXTRACT_DB(ch_ssuextract)

    INSTALL_RFAM_DB(ch_infernal)

    emit:
    checkm1    = INSTALL_CHECKM1_DB.out.db
    checkm2    = INSTALL_CHECKM2_DB.out.db
    gtdbtk     = INSTALL_GTDBTK_DB.out.db
    busco      = INSTALL_BUSCO_DB.out.db
    gvclass    = INSTALL_GVCLASS_DB.out.db
    ssuextract = INSTALL_SSUEXTRACT_DB.out.db
    rfam       = INSTALL_RFAM_DB.out.db
}


// ============================================================================
// WORKFLOW: INSTALL_ALL
// ============================================================================


workflow INSTALL_ALL {

    main:

    // Initialize all Conda environments concurrently.
    INSTALL_ENVS()

    // Each completed environment immediately emits its name.
    // Database installers receive their individual ready signals.
    // in most cases this is the conda environment installation,
    // for GVClass and SSUextract it is also the software directory
    // installation.
    ch_ready = INSTALL_ENVS.out.ready
        .mix(INSTALL_ENVS.out.software_ready)

    INSTALL_DATABASES(ch_ready)

    emit:
    environments = INSTALL_ENVS.out.ready
    checkm1      = INSTALL_DATABASES.out.checkm1
    checkm2      = INSTALL_DATABASES.out.checkm2
    gtdbtk       = INSTALL_DATABASES.out.gtdbtk
    busco        = INSTALL_DATABASES.out.busco
    gvclass      = INSTALL_DATABASES.out.gvclass
    ssuextract   = INSTALL_DATABASES.out.ssuextract
    rfam         = INSTALL_DATABASES.out.rfam
}

