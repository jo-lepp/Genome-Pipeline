process BUSCO {
    tag "$root"
    label 'busco'
    conda "bioconda::busco"
    publishDir { "${params.outdir}/${root}_busco" }, mode: params.publish_mode

    input:
    tuple val(root), path(derep_dir)

    output:
    tuple val(root), path("busco_out"),                      emit: busco_dir
    tuple val(root), path("busco_out/batch_summary.txt"),    emit: summary

    script:
    """
    busco -i "${derep_dir}" -o busco_out -l ${params.busco_lineage} -m genome \\
        --cpu ${params.threads_busco} --skip_bbtools
    """
}

workflow BUSCO_QC {
    take:
    derep_dir

    main:
    BUSCO(derep_dir)

    emit:
    busco_dir = BUSCO.out.busco_dir
    summary   = BUSCO.out.summary
}
