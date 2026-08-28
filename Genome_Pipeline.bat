Bootstrap: docker
From: continuumio/miniconda3
%post -c /bin/bash
    # Update and install system wide software here
    apt-get update

    # installing dependencies 
    conda install python=3.10
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # install conda stuff here
    conda install flye
    conda install minimap2
    conda install bioconda::samtools

    conda install bioconda::bedtools
    conda install bbmap
    conda install -c bioconda metabat2
    conda install -c conda-forge boost-cpp=1.85.0
    conda install -c bioconda myloasm mylotools
    conda install -c conda-forge -c bioconda semibin
    conda install bioconda::maxbin2
    conda install bioconda::concoct
    conda install -c bioconda -c conda-forge gsl=2.5
    conda install -c conda-forge pigz
    conda install -c bioconda miniprot
    pip install remag
    conda install -c conda-forge -c bioconda semibin
    conda install flowcraft #maxbin2