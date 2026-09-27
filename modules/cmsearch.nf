
// ---- Rfam RF01960 (SSU_rRNA_eukarya) covariance model, downloaded once ----
process CM_MODEL {
    label 'small'
    storeDir "${params.db_dir}/rfam_cm"

    output:
    path "RF01960.cm", emit: cm_file

    script:
    """
    wget -O RF01960.cm "${params.cm_model_url}"
    """
}

// ---- cmsearch + extraction, once per sample ----
process CMSEARCH {
    tag "$root"
    label 'cmsearch'
    conda "bioconda::infernal bioconda::samtools"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly)
    path cm_file

    output:
    tuple val(root), path("${root}_extracted_ssu.fa"), emit: ssu_fa
    path "${root}_ssu_hits.tblout", emit: tblout

    script:
    """
    cmsearch --cpu ${params.threads_cmsearch} --tblout ${root}_ssu_hits.tblout "${cm_file}" "${mylo_assembly}"

    # NOTE: the original fast_mode.def also runs `esl-sfetch --index` here,
    # but the actual extraction below uses `samtools faidx` (which builds
    # its own .fai index automatically) - esl-sfetch's index is never
    # actually read. Kept here only for exact parity with the source
    # script; removing it would not change the output.
    esl-sfetch --index "${mylo_assembly}"

    grep -v "^#" ${root}_ssu_hits.tblout \\
        | awk '{print \$1":"\$8"-"\$9}' \\
        | xargs samtools faidx "${mylo_assembly}" > ${root}_extracted_ssu.fa
    """
}

/*
 * Subworkflow: downloads the CM once, runs cmsearch + extraction per sample.
 */
workflow CMSEARCH_EUK {
    take:
    assembly   // tuple(root, mylo_assembly) - e.g. MYLOASM.out.assembly

    main:
    CM_MODEL()

    // .first(): same reasoning as every other one-time DB download in this
    // pipeline - built once, reused by every CMSEARCH call.
    CMSEARCH(assembly, CM_MODEL.out.cm_file.first())

    emit:
    ssu_fa = CMSEARCH.out.ssu_fa
    tblout = CMSEARCH.out.tblout
}