
// ---- one-time clone + pixi env build + reference DB setup ----
process GVCLASS_SETUP {
    label 'gvclass'
    conda "conda-forge::git conda-forge::pixi"
    storeDir "${params.db_dir}"

    output:
    path "gvclass_install", emit: install_dir

    script:
    """
    git clone ${params.gvclass_repo_url} gvclass_install
    cd gvclass_install
    pixi install --frozen
    pixi run setup-db
    """
}

// ---- .fa -> .fna copy/rename, same shape as the cp loops before each pixi run gvclass call ----
process FNA_PREP {
    tag "${root}:${label}"
    label 'small'
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(src_dir), val(pub_subpath)

    output:
    tuple val(root), val(label), path("fna"), emit: fna_dir
    path "skipped_bins_${root}_${label}.txt", emit: skipped, optional: true

    script:
    def max_bytes = params.gvclass_max_bin_size_mb * 1024 * 1024
    def skip_file = "skipped_bins_${root}_${label}.txt"
    """
    mkdir -p fna
    : > ${skip_file}

    for f in "${src_dir}"/*.fa; do
        [ -e "\$f" ] || continue

        size=\$(stat -c%s "\$f" 2>/dev/null || stat -f%z "\$f")

        if [ "\$size" -gt ${max_bytes} ]; then
            printf "%s\\t%s\\t%s\\t%s bytes\\texceeds ${params.gvclass_max_bin_size_mb}MB limit\\n" \\
                "${root}" "${label}" "\$(basename "\$f")" "\$size" >> ${skip_file}
            continue
        fi

        cp "\$f" "fna/\$(basename "\${f%.fa}.fna")"
    done

    if [ ! -s ${skip_file} ]; then
        rm -f ${skip_file}
    fi
    """
}

// ---- gather every bin set's skip report into one pipeline-wide summary ----
process COLLECT_SKIPPED {
    label 'small'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path skipped_files

    output:
    path "gvclass_skipped_bins_summary.tsv"

    script:
    """
    printf "root\\tlabel\\tbin\\tsize\\treason\\n" > gvclass_skipped_bins_summary.tsv
    cat ${skipped_files} >> gvclass_skipped_bins_summary.tsv 2>/dev/null || true
    """
}

// ---- pixi run gvclass, once per bin set (metabat/remag/self/soil/global) ----
process GVCLASS {
    tag "${root}:${label}"
    label 'gvclass'
    conda "conda-forge::pixi"
    publishDir { "${params.outdir}/${root}_gvclass" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(fna_dir)
    path install_dir

    output:
    tuple val(root), val(label), path("${label}"), emit: gvclass_dir

    script:
    """
    cd "${install_dir}"
    pixi run gvclass "\$OLDPWD/${fna_dir}" -o "\$OLDPWD/${label}" -t ${params.threads_gvclass} --tree-method fasttree
    """
}

/*
 * Subworkflow: builds GVClass once, preps .fna inputs for all five bin
 * sets, then classifies each. Unlike dRep/CheckM2/GTDB-Tk, GVClass DOES
 * run on REMAG bins (same asymmetry as BUSCO) - metabat uses its
 * dereplicated genomes, remag uses its raw bins, and self/soil/global
 * reuse the SAME unzipped bins BUSCO_QC already produced.
 */
workflow GVCLASS_TAXONOMY {
    take:
    derep_dir
    remag_bins_dir
    unzipped_dir

    main:
    GVCLASS_SETUP()

    ch_derep_metabat = derep_dir.filter { root, label, dir -> label == 'metabat' }

    ch_fna_metabat = ch_derep_metabat.map { root, label, dir -> tuple(root, label, dir, "${root}_metabat/fna") }
    ch_fna_remag   = remag_bins_dir.map   { root, dir       -> tuple(root, 'remag', dir, "${root}_remag/fna") }
    ch_fna_semibin = unzipped_dir.map     { root, label, dir -> tuple(root, label, dir, "${root}_semibin/${label}_model/fna") }

    ch_fna_in = ch_fna_metabat.mix(ch_fna_remag).mix(ch_fna_semibin)
    FNA_PREP(ch_fna_in)

    GVCLASS(FNA_PREP.out.fna_dir, GVCLASS_SETUP.out.install_dir.first())

    // gather whatever skip reports exist (possibly none) into one summary
    COLLECT_SKIPPED(FNA_PREP.out.skipped.collect())

    emit:
    gvclass_dir = GVCLASS.out.gvclass_dir
}