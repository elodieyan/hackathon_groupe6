
## y a deux manières d'installer les packages, perso c'est la qui marche  

##manière 1
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("DESeq2")
  

##manière 2
#install.packages("htmltools")
#library(htmltools)
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

##puis faut faire ça
#library( "DESeq2")
#library(ggplot2)

#browseVignettes("DESeq2") #pour accéder à la documentation

##dans la docs ils disent de plutôt faire un read.csv ou read.delim mais moi ça marche pas ils me disent "plus de colonnes que de noms de colonnes"
##donc j'avais fait un read.table qui lui marche bien 
file = read.delim("C:/Users/Administrator/Desktop/AMI2B/Hackathon/SRR628582.counts",sep="\t",  header=T, dec=",")
#file = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/SRR628582.counts",sep="\t",  header=T, dec=",")


##je ne sais pas quoi mettre dans colData =, j'ai une erreur à chaque fois 
dds <- DESeqDataSetFromMatrix(countData = file,
                              colData = Geneid,
                              design= ~condition)
dds <- DESeq(dds)
file <- results(dds)


