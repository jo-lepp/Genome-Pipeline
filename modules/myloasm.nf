process MYLOASM {
    tag "$root"
    label 'assembly'
    conda "bioconda::myloasm"
    publishDir { "${params.outdir}/${root}_myloasm" }, mode: params.publish_mode
    
    input:
    tuple val(root), path(input_fgz)

    output:
    tuple val(root), path("assembly_primary.fa"), emit: assembly
    path "*", emit: mylo_bins_dir

    script:
    def preset_flag = params.read_type == 'hifi' ? '--hifi' : '--nano-r10'
    """
    myloasm "${input_fgz}" -o . -t ${params.threads_myloasm} ${preset_flag}

    echo "Myloasm assembly done (${params.read_type})."
    """
}