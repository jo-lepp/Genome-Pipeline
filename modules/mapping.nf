process MAPPING {
    tag "$root"
    label 'mapping'
    conda "bioconda::minimap2 bioconda::samtools"
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: "*.{sam,bam,bai}"

    input:
    tuple val(root), path(mylo_assembly), path(input_fgz)

    output:
    tuple val(root), path("${root}_sorted.bam"), path("${root}_sorted.bam.bai"), emit: bam
    path "${root}_.sam", emit: sam

    script:
    def minimap_preset = params.read_type == 'hifi' ? 'map-hifi' : 'map-ont'
    """
    minimap2 --eqx -t ${params.threads_map} -a -x ${minimap_preset} "${mylo_assembly}" "${input_fgz}" > "${root}_.sam"
    samtools sort -@${params.threads_map} -o "${root}_sorted.bam" "${root}_.sam"
    samtools index "${root}_sorted.bam"
    """
}