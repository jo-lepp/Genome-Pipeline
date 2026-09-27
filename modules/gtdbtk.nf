
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
    gtdbtk classify_wf --genome_dir "${derep_dir}" --out_dir gtdbtk -x .${ext}
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
    derep_dir   // tuple(root, label, dereplicated_genomes dir) - from DEREPLICATION.out.derep_dir

    main:
    GTDBTK_DB()

    ch_gtdbtk_in = derep_dir.map { root, label, dir ->
        def ext = (label == 'metabat') ? 'fa' : 'fa.gz'
        def sub = (label == 'metabat') ? "${root}_metabat/gtdbtk" : "${root}_semibin/${label}_model/gtdbtk"
        tuple(root, label, dir, ext, sub)
    }

    // same reasoning as CHECKM2_DB: this ran once, so .first() makes it
    // reusable across every GTDBTK invocation instead of being drained
    // after the first one.
    GTDBTK(ch_gtdbtk_in, GTDBTK_DB.out.db_dir.first())

    emit:
    gtdbtk_dir = GTDBTK.out.gtdbtk_dir
}