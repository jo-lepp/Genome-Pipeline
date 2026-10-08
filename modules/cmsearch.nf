/*
 * NEW default fast-mode tool: SSUextract (https://github.com/NeLLi-team/ssuextract)
 *
 * A full Nextflow pipeline in its own right (pixi-managed, same pattern as
 * GVClass): detects 16S/18S rRNA genes via bundled Infernal covariance
 * models (RF00177 -> 16S index, RF01960 -> 18S index), then assigns
 * taxonomy via BLAST against SILVA+PR2 (curated profile) or SILVA+PR2+IMG
 * (img profile), with an optional tree-based classification step.
 * Archaea are covered through SILVA's 16S taxonomy, not a separate model.
 */

process SSUEXTRACT_SETUP {
    label 'ssuextract'
    conda "conda-forge::git conda-forge::pixi"
    storeDir "${params.db_dir}"

    output:
    path "ssuextract_install", emit: install_dir

    script:
    """
    git clone ${params.ssuextract_repo_url} ssuextract_install
    cd ssuextract_install

    # SSUextract is its own nested Nextflow pipeline with its own internal
    # per-process time limit (config/base.config), separate from anything
    # in OUR nextflow.config - raising our own threads_cmsearch has no
    # effect on it. Strip the time limit so long BLAST/annotation steps
    # on large assemblies aren't killed mid-run by the inner pipeline's
    # own unrelated timeout.
    sed -i '/^[[:space:]]*time[[:space:]]*=/s/^/# /' config/base.config

    pixi install --frozen
    pixi run setup --database_profile ${params.ssuextract_db_profile}
    """
}

process SSU_EXTRACT {
    tag "$root"
    label 'ssuextract'
    conda "conda-forge::pixi"
    publishDir { "${params.outdir}/${root}_ssuextract" }, mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly)
    path install_dir

    output:
    tuple val(root), path("results"), emit: results_dir
    path "results/cmsearch_summary.tsv", emit: summary, optional: true
    path "results/extracted/*.fna", emit: extracted_fastas, optional: true

    script:
    def tree_flag = params.ssuextract_tree_classification ? '--tree_classification' : ''
    """
    cd "${install_dir}"
    pixi run ssuextract -q "\$OLDPWD/${mylo_assembly}" --outdir "\$OLDPWD/results" ${tree_flag}
    """
}

workflow SSU_EXTRACT_ALL {
    take:
    assembly   // tuple(root, mylo_assembly)

    main:
    SSUEXTRACT_SETUP()
    SSU_EXTRACT(assembly, SSUEXTRACT_SETUP.out.install_dir.first())

    emit:
    results_dir = SSU_EXTRACT.out.results_dir
    summary     = SSU_EXTRACT.out.summary
}

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
