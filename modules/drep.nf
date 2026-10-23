/*
 * NOTE: dRep is never run on REMAG bins in the original script - only
 * metabat2 and the three SemiBin2 models get dereplicated.
 *
 * IMPORTANT: dRep's default filtering step internally shells out to
 * checkM v1 (a separate, older tool from CheckM2, despite the similar
 * name) to score completeness/contamination before clustering. checkM1
 * needs its own reference database registered via `checkm data setRoot`
 * before it will run at all - without it, dRep fails at Step 1 with
 * "!!! checkM failed !!!" and silently produces zero dereplicated
 * genomes, even though dRep itself exits with status 0.
 */

// ---- checkM v1 reference database, downloaded once ----
process CHECKM1_DB {
    label 'binning'
    conda "bioconda::drep conda-forge::python=3.10 conda-forge::setuptools=75.8.2"
    storeDir "${params.db_dir}/checkm1_db"

    output:
    path "checkm_data", emit: db_dir

    script:
    """
    wget -O checkm_data.tar.gz "${params.checkm1_data_url}"
    mkdir -p checkm_data
    tar -xzf checkm_data.tar.gz -C checkm_data
    rm checkm_data.tar.gz
    """
}

process DREP {
    tag "${root}:${label}"
    label 'binning'
    conda "bioconda::drep conda-forge::python=3.10 conda-forge::setuptools=75.8.2"
    publishDir { "${params.outdir}/${pub_subpath}" }, mode: params.publish_mode

    input:
    tuple val(root), val(label), path(genomes, stageAs: 'genomes_in/*'), val(ext), val(pub_subpath)
    path checkm1_db

    output:
    tuple val(root), val(label), path("drep/dereplicated_genomes"), emit: derep_dir

    script:
    def quality_flag = params.drep_ignore_genome_quality ? '--ignoreGenomeQuality' : ''
    def checkm_setup = params.drep_ignore_genome_quality ? '' : "checkm data setRoot \"\$(pwd)/${checkm1_db}\""
    """
    ${checkm_setup}
    # bins may arrive gzipped (SemiBin2) or plain (metabat2/REMAG) - normalize
    # to plain fasta before dRep ever sees them, so every downstream consumer
    # (BUSCO, CheckM2, GTDB-Tk, GVClass) can assume uniform, uncompressed input
    mkdir -p plain_bins
    for f in bin_*.fa; do
        if file "\$f" | grep -q gzip; then
            gunzip -c "\$f" > "plain_bins/\$(basename "\$f")"
        else
            cp "\$f" "plain_bins/"
        fi
    done
    
    dRep dereplicate drep_out -g plain_bins/*.fa ${quality_flag}
    cp -r drep_out/dereplicated_genomes .
    """
}

/*
 * Subworkflow: takes the bin channels emitted by BINNING and runs dRep
 * on whichever of metabat2/semibin2 actually produced bins.
 */
workflow DEREPLICATION {
    take:
    metabat_bins   // tuple(root, [fa files])          - from BINNING.out.metabat_bins
    semibin_bins   // tuple(root, model, [fa.gz files]) - from BINNING.out.semibin_bins

    main:
    CHECKM1_DB()

    ch_drep_metabat = metabat_bins.map { root, fa ->
        tuple(root, 'metabat', fa, 'fa', "${root}_metabat/drep")
    }
    ch_drep_semibin = semibin_bins.map { root, model, fagz ->
        tuple(root, model, fagz, 'fa.gz', "${root}_semibin/${model}_model/drep")
    }

    ch_drep_in = ch_drep_metabat.mix(ch_drep_semibin)
    DREP(ch_drep_in, CHECKM1_DB.out.db_dir.first())

    emit:
    derep_dir = DREP.out.derep_dir
}
