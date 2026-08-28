Bootstrap: docker
From: continuumio/miniconda3
%post -c /bin/bash
    # Update and install system wide software here
    apt-get update

    # installing dependencies
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # install conda stuff here
    conda create --name assemblers flye myloasm mylotools 
    conda create --name tools minimap2 samtools bedtools bbmap pigz
    conda create --name binners metabat2 semibin concoct python=3.10 boost-cpp=1.85.0

%arguments
    #have variables in the environment section instead?
    input_bam="/home/user/data"
    input_adapters="/home/user/data"
    output_path="/home/user/data"
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
    conda activate tools
    prefix="${input_bam%.bam}"
    root="${prefix##*/}"
    out_bfq="${output_path}/${root}.bam.fastq"
    bedtools bamtofastq -i "$input_bam" -fq "$out_bfq"
    pigz "$out_bfq"

    #Continue to filter

    #have ktriptips and k and edist be inputs?
    #trimming adapters
    out_trim_fgz="${output_path}/${root}.fastq.gz"
    #change this to filtered fgz (by this i mean out_fgz as input below)
    bbduk.sh in="${out_bfq}.gz" out="$out_trim_fgz" ktrimtips=50 ref="$input_adapters" k=21 edist=1

    echo "Processing done."
    conda deactivate
#ASSEMBLING WITH FLYE
    conda activate assemblers
    flye_bins="${output_path}/${root}_flye_bins/"
    flye --meta -t 10 --pacbio-hifi "$out_trim_fgz" -o "$flye_bins"
    flye_assembly="${flye_bins}/assembly.fasta"
        #flye makes assembly.fasta

    echo "Flye assembly done."
    conda deactivate
    #minimap, samtools, metabat2 for flye
    conda activate tools
    sam_filepath="${output_path}/${root}_fminimap.sam"
    minimap2 --eqx -t 32 -a -x map-hifi "$flye_assembly" "$out_trim_fgz" > "$sam_filepath"
    bam_path="${output_path}/${root}_fsort.bam"
    samtools view -b -F 0x800 -F 0x100 -T "$flye_assembly" "$sam_filepath" | samtools sort -@ 32 -o "$bam_path"
    samtools index "$bam_path"
    depth_path="${output_path}/${root}_flye_depth_test.txt"
    jgi_summarize_bam_contig_depths  --outputDepth "$depth_path" "$bam_path"
    conda deactivate
    conda activate binners
    metabat_path="${output_path}/${root}_flye_metabat_bins"
    metabat2 -i "$flye_assembly" -a "$depth_path" -o "$metabat_path"

    echo "Flye metabat2 binning done."
    
    conda create concoct_env python=3 concoct
    conda activate concoct_env
    #concoct, flye
    concoct_output="${output_path}/${root}_flye_concoct"
    concoct_bed="${concoct_output}/${root}_concoct.bed"
    concoct_fa="${concoct_output}/${root}_concoct.fa"
    cut_up_fasta.py "${flye_assembly}" -c 10000 -o 0 --merge_last -b "${concoct_bed}" > "${concoct_fa}"
    concoct_table="${concoct_output}/${root}_concoct_table.tsv"
    concoct_coverage_table.py "${concoct_bed}" "$bam_path" > "${concoct_table}"
    concoct --composition_file "${concoct_fa}" --coverage_file "${concoct_table}" -b "${concoct_output}"
    concoct_csv="${concoct_output}/${root}_clustering.csv"
    concoct_merged="${concoct_output}/${root}_concoct_merged.csv"
    merge_cutup_clustering.py "${concoct_csv}" > "${concoct_merged}"
    concoct_fasta_dir="${concoct_output}/${root}_fasta_bins"
    mkdir "${concoct_fasta_dir}"
    extract_fasta_bins.py "$flye_assembly" "${concoct_merged}" --output_path "${concoct_fasta_dir}"

    echo "Flye concoct binning done."

    conda deactivate
    

#ASSEMBLING WITH MYLOASM
    conda activate assemblers
    mylo_bins="${output_path}/${root}_mylo_bins"
    myloasm "$out_trim_fgz" -o "$mylo_bins" -t 10 --hifi
    mylo_assembly="${mylo_bins}/${root}_assembly_primary.fa"
        #mylo makes primary_assembly.fa
    conda deactivate

#minimap, samtools, metabat2 for myloasm
    conda activate tools
    sam_filepath="${output_path}/${root}_mylominimap.sam"
    minimap2 --eqx -t 32 -a -x map-hifi "$mylo_assembly" "$out_trim_fgz" > "$sam_filepath"
    bam_path="${output_path}/${root}_mysort.bam"
    samtools view -b -F 0x800 -F 0x100 -T "$mylo_assembly" "$sam_filepath" | samtools sort -@ 32 -o "$bam_path"
    samtools index "$bam_path"
    depth_path="${output_path}/${root}_mylo_depth_test.txt"
    jgi_summarize_bam_contig_depths  --outputDepth "$depth_path" "$bam_path"
    conda deactivate
    conda activate binners
    metabat_path="${output_path}/${root}_mylo_metabat_bins"
    metabat2 -i "$mylo_assembly" -a "$depth_path" -o "$metabat_path"

    echo "Myloasm assembly done."
    