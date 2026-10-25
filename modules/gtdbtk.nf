
// ---- reference package download, run once ----
process GTDBTK_DB {
    label 'taxonomy'
    conda "bioconda::gtdbtk"
    storeDir "${params.db_dir}/gtdbtk_db"

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

// ---- gtdbtk classify_wf, run once per sample on dRep's combined output ----
process GTDBTK {
    tag "$root"
    label 'taxonomy'
    conda "bioconda::gtdbtk"
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
    gtdbtk classify_wf --genome_dir "${derep_dir}" --out_dir gtdbtk -x fa \\
        --cpus ${params.threads_gtdbtk} \\
        --pplacer_cpus 1 \\
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
    GTDBTK_DB()
    GTDBTK(derep_dir, GTDBTK_DB.out.db_dir)

    emit:
    classify_dir = GTDBTK.out.classify_dir
    summary      = GTDBTK.out.summary
}
