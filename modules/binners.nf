process DEPTH {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"

    input:
    tuple val(root), path(bam), path(bai)

    output:
    tuple val(root), path("${root}_depth.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth "${root}_depth.txt" "${bam}"
    """
}

process METABAT2 {
    tag "$root"
    label 'binning'
    conda "bioconda::metabat2"
    publishDir { "${params.outdir}/${root}_metabat" }, mode: params.publish_mode

    input:
    tuple val(root), path(assembly), path(depth_path)

    output:
    tuple val(root), path("bins/*.fa"), emit: bins

    script:
    """
    mkdir -p bins
    metabat2 -i "${assembly}" -a "${depth_path}" -o bins/metabat
    """
}

process SEMIBIN2 {
    tag "${root}:${model}"
    label 'binning'
    conda "bioconda::semibin"
    publishDir { "${params.outdir}/${root}_semibin/${model}_model" }, mode: params.publish_mode

    input:
    tuple val(root), path(assembly), path(bam), path(bai), val(model), val(mode_flag)

    output:
    tuple val(root), path("output_bins/*.fa.gz"), emit: bins

    script:
    """
    SemiBin2 single_easy_bin --sequencing-type long_read ${mode_flag} -i "${assembly}" -b "${bam}" -o .
    """
}
