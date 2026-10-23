
// ---- reference package download, run once ----
process GTDBTK_DB {
    label 'taxonomy'
    conda "bioconda::gtdbtk"
    storeDir "${params.db_dir}/gtdbtk_db"   // <-- replaces publishDir

    output:
    path "gtdbtk_data", emit: db_dir

    script:
    """
    mkdir -p gtdbtk_data
    wget -O gtdbtk_data.tar.gz "${params.gtdbtk_data_url}"
    tar -xzf gtdbtk_data.tar.gz -C gtdbtk_data --strip-components=1
    rm gtdbtk_data.tar.gz
    """
}

// ---- gtdbtk classify_wf, run once per bin set (metabat/self/soil/global) ----
process GTDBTK {
    tag "${root}:${label}"
    label 'taxonomy'
    conda "bioconda::gtdbtk"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(derep_dir), val(ext), val(pub_subpath)
    path db_dir

    output:
    tuple val(root), val(label), path("gtdbtk"), emit: gtdbtk_dir

    script:
    """
    export GTDBTK_DATA_PATH="\$(pwd)/${db_dir}"
    mkdir -p pplacer_scratch
    gtdbtk classify_wf --genome_dir "${derep_dir}" --out_dir gtdbtk -x .${ext} \\
        --cpus ${params.threads_gtdbtk} \\
        --pplacer_cpus 1 \\
        --scratch_dir pplacer_scratch
    echo "bacterial taxonomy done (${label})."
    """
}

/*
 * Subworkflow: downloads the reference package once, then classifies
 * whatever DEREPLICATION handed back (metabat + whichever SemiBin2
 * models ran) - same as CHECKM2_QC, nothing extra to gate here since
 * a disabled binner already means an empty upstream channel.
 */
workflow GTDBTK_CLASSIFY {
    take:
    derep_dir   // tuple(root, dir) - was tuple(root, label, dir)

    main:
    GTDBTK_DB()
    GTDBTK(derep_dir, GTDBTK_DB.out.db_dir.first())

    emit:
    classify_dir = GTDBTK.out.classify_dir
    summary      = GTDBTK.out.summary
}
