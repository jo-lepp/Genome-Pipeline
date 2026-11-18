// ---- gtdbtk classify_wf, run once per sample on dRep's combined output ----
process GTDBTK {
    tag "$root"
    label 'taxonomy'
    conda "${projectDir}/envs/gtdbtk.yml"
    publishDir { "${params.outdir}/${root}_gtdbtk" }, mode: params.publish_mode

    input:
    tuple val(root), path(derep_dir)
    path db_dir

    output:
    tuple val(root), path("gtdbtk"),               emit: classify_dir
    tuple val(root), path("gtdbtk_summary.tsv"),   emit: summary

    script:
    """
    export GTDBTK_DATA_PATH="\$(pwd)/${db_dir}"
    mkdir -p pplacer_scratch
    # GTDB gets all of the cpus provided to the task; 
    # pplacer receives half of them. 
    gtdbtk classify_wf --genome_dir "${derep_dir}" --out_dir gtdbtk -x fa \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${Math.max(1, task.cpus.intdiv(2))} \\
        --scratch_dir pplacer_scratch

    # classify_wf writes bac120/ar53 summaries separately, whichever apply -
    # combine into one file so downstream (REPORT) only has one thing to read
    : > gtdbtk_summary.tsv
    header_written=false
    for f in gtdbtk/gtdbtk.bac120.summary.tsv gtdbtk/gtdbtk.ar53.summary.tsv; do
        if [ -f "\$f" ]; then
            if [ "\$header_written" = false ]; then
                cp "\$f" gtdbtk_summary.tsv
                header_written=true
            else
                tail -n +2 "\$f" >> gtdbtk_summary.tsv
            fi
        fi
    done

    echo "bacterial/archaeal taxonomy done."
    """
}

workflow GTDBTK_CLASSIFY {
    take:
    derep_dir   // tuple(root, dir)

    main:
    ch_db = Channel.value(file("${params.db_dir}/gtdbtk_db/gtdbtk_data", checkIfExists: true))
    GTDBTK(derep_dir, ch_db)

    emit:
    classify_dir = GTDBTK.out.classify_dir
    summary      = GTDBTK.out.summary
}
