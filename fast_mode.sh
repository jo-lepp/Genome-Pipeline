Bootstrap: localimage
From: /home/jalepper/genome_pipeline/prototyping/take_two/env.sif

%arguments
    #have variables in the environment section instead?
    input_fgz="/home/user/data/input.fastq.gz"
    output_path="/home/user/data/output"
    root="sample"

%environment
    export PATH="/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:$PATH"
    export CONDA_DEFAULT_ENV=base

%runscript
    #!/bin/bash
    set -eo pipefail
    . /opt/conda/etc/profile.d/conda.sh
    conda activate base
#PROCESSING
    prefix="${input_fgz%.fastq.gz}"
    root="${prefix##*/}"
    
#ASSEMBLING WITH MYLOASM
    conda activate assemblers
    mylo_bins="${output_path}/${root}_myloasm"
    myloasm "$input_fgz" -o "$mylo_bins" -t 10 --hifi
    mylo_assembly="${mylo_bins}/assembly_primary.fa"
    conda deactivate

#minimap, samtools
    conda activate tools
    sam_filepath="${output_path}/${root}_.sam"
    minimap2 --eqx -t 32 -a -x map-hifi "$mylo_assembly" "$input_fgz" > "$sam_filepath"
    bam_path="${output_path}/${root}_sorted.bam"
    samtools sort -@ 32 -o "$bam_path" "$sam_filepath"
    samtools index "$bam_path"
    conda deactivate
    echo "Myloasm assembly done."

#binning (metabat2, semibin, remag)
    conda activate binners
    depth_path="${output_path}/${root}_depth.txt"
    jgi_summarize_bam_contig_depths  --outputDepth "$depth_path" "$bam_path"

    mkdir "${output_path}/${root}_metabat"
    mkdir "${output_path}/${root}_metabat/bins"
    metabat_path="${output_path}/${root}_metabat/bins"
    metabat_output="${output_path}/${root}_metabat/bins/metabat"
    metabat2 -i "$mylo_assembly" -a "$depth_path" -o "$metabat_output" --unbinned
    mv ${metabat_path}/*.txt ${output_path}/${root}_metabat
    echo "metabat2 binning done"

    mkdir "${output_path}/${root}_semibin"
    mkdir "${output_path}/${root}_semibin/soil_model"
    soil_path="${output_path}/${root}_semibin/soil_model"
    mkdir "${output_path}/${root}_semibin/global_model"
    global_path="${output_path}/${root}_semibin/global_model"

    SemiBin2 single_easy_bin --sequencing-type long_read --environment soil -i $mylo_assembly -b $bam_path -o $soil_path
    SemiBin2 single_easy_bin --sequencing-type long_read --environment global -i $mylo_assembly -b $bam_path -o $global_path
    echo "Semibin2 binning done"

    conda deactivate

#DEREPLICATION
    conda activate drep
    dRep dereplicate "${output_path}/${root}_metabat/drep" -g ${metabat_path}/*.fa
    dRep dereplicate "${soil_path}/drep" -g ${soil_path}/output_bins/*.fa.gz
    dRep dereplicate "${global_path}/drep" -g ${global_path}/output_bins/*.fa.gz
    echo "dereplication done"
    conda deactivate

#BACTERIAL QC
    conda activate checkm2
    mkdir ${output_path}/checkm2_db
    db_path="${output_path}/checkm2_db"
    checkm2 database --download --path $db_path --no_write_json_db
    checkm2 predict --threads 10 --input "${output_path}/${root}_metabat/drep/dereplicated_genomes" --extension .fa --output_directory "${output_path}/${root}_metabat/checkm2" --database_path ${db_path}/CheckM2_database/uniref100.KO.1.dmnd
    checkm2 predict --threads 10 --input "${soil_path}/drep/dereplicated_genomes" --extension .fa.gz --output_directory "${soil_path}/checkm2" --database_path ${db_path}/CheckM2_database/uniref100.KO.1.dmnd
    checkm2 predict --threads 10 --input "${global_path}/drep/dereplicated_genomes" --extension .fa.gz --output_directory "${global_path}/checkm2" --database_path ${db_path}/CheckM2_database/uniref100.KO.1.dmnd
    echo "checkm2 done"
    conda deactivate

#EUKARYOTIC QC
    conda activate busco

    soil_unzipped="${soil_path}/output_bins_unzipped"
    mkdir -p "$soil_unzipped"
    for f in "${soil_path}/drep/dereplicated_genomes/"*.fa.gz; do
        gunzip -c "$f" > "${soil_unzipped}/$(basename "${f%.gz}")"
    done
    global_unzipped="${global_path}/output_bins_unzipped"
    mkdir -p "$global_unzipped"
    for f in "${global_path}/drep/dereplicated_genomes/"*.fa.gz; do
        gunzip -c "$f" > "${global_unzipped}/$(basename "${f%.gz}")"
    done
    busco -i ${output_path}/{root}_metabat/drep/dereplicated_genomes --out_path "${output_path}/${root}_metabat" -o busco -m genome -f -l eukaryota_odb10
    busco -i "${soil_path}/output_bins_unzipped" --out_path ${soil_path} -o busco -m genome -f -l eukaryota_odb10
    busco -i "${global_path}/output_bins_unzipped" --out_path $global_path -o busco -m genome -f -l eukaryota_odb10
    echo "busco done."
    conda deactivate

#CMSEARCH

    
    echo "Fast_mode pipeline done at [$(date '+%H:%M:%S')]."

    
