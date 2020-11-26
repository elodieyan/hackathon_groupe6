BootStrap: docker
From: r-base:4.0.2

%labels
Maintainer Guillaume Studer
base.image="r-base:4.0.2"
version="1"
software="EnhencedVolcano+DESeq2+FactoMineR+factoextra"
%help
Please faites que ca fonctionne
%post
apt-get update
apt-get install -y procps libssl-dev libcurl4-gnutls-dev curl git libopenmpi-dev openmpi-bin openmpi-doc libxml2-dev
R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}'
R -e 'BiocManager::install("DESeq")'
R -e 'install.packages("FactoMineR")'
R -e 'install.packages("factoextra")'
