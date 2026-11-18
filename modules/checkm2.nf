process CHECKM2 {
    tag "$root"
    label 'checkm2'
    conda "${projectDir}/envs/checkm2.yml"
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
    checkm2 predict --threads ${task.cpus} --input "${derep_dir}" -x fa \\
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
    ch_db = Channel.value(file("${params.db_dir}/checkm2_data/checkm2_db", checkIfExists: true))
    CHECKM2(derep_dir, ch_db)

    emit:
    checkm2_dir     = CHECKM2.out.checkm2_dir
    quality_report  = CHECKM2.out.quality_report   // feeds REPORT directly
}
