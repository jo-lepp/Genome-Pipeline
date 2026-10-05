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
// ---- .fa -> .fna copy/rename, with oversized bins subsampled via `gt splitfasta` ----
// Previously: bins over params.gvclass_max_bin_size_mb were excluded entirely.
// Now: they're split into multiple smaller files (each under the size limit)
// using `gt splitfasta -targetsize`, and ALL resulting parts are still sent
// to GVClass - logged as "subsampled", not dropped.
process FNA_PREP {
    tag "${root}:${label}"
    label 'small'
    conda "bioconda::genometools-genometools"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(src_dir), val(pub_subpath)

    output:
    tuple val(root), val(label), path("fna"), emit: fna_dir
    path "subsampled_bins_${root}_${label}.txt", emit: subsampled, optional: true

    script:
    def max_bytes = params.gvclass_max_bin_size_mb * 1024 * 1024
    def report = "subsampled_bins_${root}_${label}.txt"
    """
    mkdir -p fna
    : > ${report}

    for f in "${src_dir}"/*.fa; do
        [ -e "\$f" ] || continue
        base=\$(basename "\$f" .fa)
        size=\$(stat -c%s "\$f" 2>/dev/null || stat -f%z "\$f")

        if [ "\$size" -le ${max_bytes} ]; then
            cp "\$f" "fna/\${base}.fna"
            continue
        fi

        size_mb=\$(( (size + 1048575) / 1048576 ))

        # Isolate this one bin in its own directory so gt splitfasta's
        # output files (whatever it names them) can be identified simply
        # as "everything in this directory" - no other files to confuse with.
        mkdir -p "split_tmp/\${base}"
        cp "\$f" "split_tmp/\${base}/\${base}.fa"
        ( cd "split_tmp/\${base}" && gt splitfasta -targetsize ${params.gvclass_max_bin_size_mb} "\${base}.fa" )
        rm -f "split_tmp/\${base}/\${base}.fa"

        part_count=0
        oversized_part_warning=""
        for part in split_tmp/\${base}/*; do
            [ -e "\$part" ] || continue
            part_count=\$((part_count + 1))
            part_size=\$(stat -c%s "\$part" 2>/dev/null || stat -f%z "\$part")
            if [ "\$part_size" -gt ${max_bytes} ]; then
                oversized_part_warning=" [WARNING: part \$part_count still exceeds limit - likely a single oversized contig that can't be split further]"
            fi
            cp "\$part" "fna/\${base}_part\${part_count}.fna"
        done

        printf "%s\\t%s\\t%s\\t%s bytes (~%s MB)\\tsubsampled into %s part(s) via gt splitfasta -targetsize ${params.gvclass_max_bin_size_mb}%s\\n" \\
            "${root}" "${label}" "\${base}.fa" "\$size" "\$size_mb" "\$part_count" "\$oversized_part_warning" >> ${report}
    done

    if [ ! -s ${report} ]; then
        rm -f ${report}
    fi
    """
}

// ---- gather every bin set's subsample report into one pipeline-wide summary ----
process COLLECT_SUBSAMPLE_REPORT {
    label 'small'
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    path subsampled_files

    output:
    path "gvclass_subsampled_bins_summary.tsv"

    script:
    """
    printf "root\\tlabel\\tbin\\tsize\\tnote\\n" > gvclass_subsampled_bins_summary.tsv
    cat ${subsampled_files} >> gvclass_subsampled_bins_summary.tsv 2>/dev/null || true
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

    COLLECT_SUBSAMPLE_REPORT(FNA_PREP.out.subsampled.collect())

    emit:
    gvclass_dir = GVCLASS.out.gvclass_dir
}