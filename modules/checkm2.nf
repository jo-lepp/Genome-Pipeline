process CHECKM2_DB {
    label 'small'
    storeDir "${params.db_dir}/checkm2_data"

    output:
    path "checkm2_db", emit: db_path

    script:
    """
    checkm2 database --download --path checkm2_db
    """
}

process CHECKM2 {
    tag "$root"
    label 'checkm2'
    conda "bioconda::checkm2"
    publishDir { "${params.outdir}/${root}_checkm2" }, mode: params.publish_mode

    input:
    tuple val(root), path(derep_dir)
    path checkm2_db

    output:
    tuple val(root), path("quality_report.tsv"), emit: quality_report
    path "checkm2_out", emit: checkm2_dir

    script:
    """
    export CHECKM2DB="\$(pwd)/${checkm2_db}/CheckM2_database/uniref100.KO.1.dmnd"
    checkm2 predict --threads ${params.threads_checkm2} --input "${derep_dir}" -x fa \\
        --output-directory checkm2_out
    cp checkm2_out/quality_report.tsv .
    """
}

/*
 * Simplified: one dRep output, one CheckM2 run per sample. No more
 * per-label iteration (metabat/soil/self/global) and no more excluding
 * remag - with only one binner active per run, there's nothing left to
 * iterate over.
 */
workflow CHECKM2_QC {
    take:
    derep_dir   // tuple(root, dir)

    main:
    CHECKM2_DB()
    CHECKM2(derep_dir, CHECKM2_DB.out.db_path.first())

    emit:
    checkm2_dir     = CHECKM2.out.checkm2_dir
    quality_report  = CHECKM2.out.quality_report   // feeds REPORT directly
}

