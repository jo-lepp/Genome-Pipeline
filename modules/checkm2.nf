// ---- database download, run once and reused by every predict call ----
process CHECKM2_DB {
    label 'qc'
    conda "bioconda::checkm2"
    storeDir "${params.db_dir}/checkm2_db"   // <-- replaces publishDir

    output:
    path "CheckM2_database/uniref100.KO.1.dmnd", emit: db

    script:
    """
    checkm2 database --download --path . --no_write_json_db
    """
}

// ---- checkm2 predict, run once per bin set (metabat/self/soil/global) ----
process CHECKM2 {
    tag "${root}:${label}"
    label 'qc'
    conda "bioconda::checkm2"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(derep_dir), val(ext), val(pub_subpath)
    path db

    output:
    tuple val(root), val(label), path("checkm2"), emit: checkm2_dir

    script:
    """
    checkm2 predict --threads ${params.threads_checkm2} \\
        --input "${derep_dir}" --extension .${ext} \\
        --output_directory checkm2 \\
        --database_path "${db}"

    echo "checkm2 done (${label})"
    """
}

/*
 * Subworkflow: downloads the DB once, then runs checkm2 predict on
 * whatever DEREPLICATION handed back (metabat + whichever SemiBin2
 * models ran) - nothing extra to gate here, since an empty upstream
 * channel (a binner that was turned off) just means nothing to predict on.
 */
workflow CHECKM2_QC {
    take:
    derep_dir   // tuple(root, label, dereplicated_genomes dir) - from DEREPLICATION.out.derep_dir

    main:
    CHECKM2_DB()

    ch_checkm2_in = derep_dir.map { root, label, dir ->
        def ext = (label == 'metabat') ? 'fa' : 'fa.gz'
        def sub = (label == 'metabat') ? "${root}_metabat/checkm2" : "${root}_semibin/${label}_model/checkm2"
        tuple(root, label, dir, ext, sub)
    }

    // .first() turns the single CHECKM2_DB output into a reusable "value"
    // channel so every CHECKM2 call below can consume it, instead of it
    // being drained after the first invocation like a normal queue channel.
    CHECKM2(ch_checkm2_in, CHECKM2_DB.out.db.first())

    emit:
    checkm2_dir = CHECKM2.out.checkm2_dir
}

