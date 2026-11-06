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
    sed -i '/^[[:space:]]*time[[:space:]]*=/d' config/base.config

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
    tuple val(root), path("results/cmsearch_summary.tsv"), emit: summary, optional: true
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
    wget -O RF01960.cm "${params.rfam_eukarya_cm_url}"
    """
}

// ---- cmsearch + extraction, once per sample ----
process CMSEARCH {
    tag "$root"
    label 'cmsearch'
    conda "bioconda::infernal bioconda::samtools bioconda::seqkit"
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    tuple val(root), path(mylo_assembly)
    path cm_file

    output:
    tuple val(root), path("${root}_extracted_ssu.fa"),         emit: ssu_fa
    tuple val(root), path("${root}_cmsearch_summary.tsv"),     emit: summary
    path "${root}_ssu_hits.tblout",                            emit: tblout

    script:
    """
    cmsearch --cpu ${params.threads_cmsearch} --tblout ${root}_ssu_hits.tblout "${cm_file}" "${mylo_assembly}"

    # cmsearch reports reverse-strand hits with seq-from > seq-to (column
    # 10 is the strand itself). samtools faidx requires start <= end and
    # silently returns an empty sequence otherwise - swap coordinates for
    # '-' strand hits before fetching, then reverse-complement the result.
    #
    # Also builds the summary.tsv REPORT expects (contig_id/taxonomy/sequence)
    # in this same loop, since the extraction already walks every hit once.
    # NOTE: a single Rfam covariance model only confirms "this is eukaryotic
    # SSU rRNA" - it has no BLAST-against-reference step like SSUextract, so
    # it genuinely cannot assign a species/taxon name. The taxonomy column
    # records that limitation honestly rather than fabricating a classification.
    : > ${root}_extracted_ssu.fa
    echo -e "contig_id\\ttaxonomy\\tsequence" > ${root}_cmsearch_summary.tsv

    grep -v "^#" ${root}_ssu_hits.tblout | awk '{print \$1, \$8, \$9, \$10}' | \\
    while read -r tgt tfrom tto strand; do
        if [ "\$strand" = "-" ]; then
            region="\${tgt}:\${tto}-\${tfrom}"
            seq=\$(samtools faidx "${mylo_assembly}" "\$region" | seqkit seq -r -p -t dna)
        else
            region="\${tgt}:\${tfrom}-\${tto}"
            seq=\$(samtools faidx "${mylo_assembly}" "\$region")
        fi
        echo "\$seq" >> ${root}_extracted_ssu.fa

        seq_only=\$(echo "\$seq" | grep -v "^>" | tr -d '\\n')
        echo -e "\${tgt}\\tEukaryota (SSU rRNA, RF01960 - no further classification)\\t\${seq_only}" >> ${root}_cmsearch_summary.tsv
    done
    """
}

workflow CMSEARCH_EUK {
    take:
    assembly   // tuple(root, mylo_assembly)

    main:
    CM_MODEL()
    CMSEARCH(assembly, CM_MODEL.out.cm_file.first())

    emit:
    ssu_fa  = CMSEARCH.out.ssu_fa
    summary = CMSEARCH.out.summary
    tblout  = CMSEARCH.out.tblout
}
