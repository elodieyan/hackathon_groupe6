



##installation des packages 

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("DESeq2")
  

##puis faut faire ça
#library( "DESeq2")
#library(ggplot2)

#browseVignettes("DESeq2") #pour accéder à la documentation

##dans la docs ils disent de plutôt faire un read.csv ou read.delim mais moi ça marche pas ils me disent "plus de colonnes que de noms de colonnes"
##donc j'avais fait un read.table qui lui marche bien 
#file = read.delim("C:/Users/Administrator/Desktop/AMI2B/Hackathon/SRR628582.counts",sep="\t",  header=T, dec=",")
file = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/SRR628582.counts",sep="\t",  header=T, dec=",")

#pour CountData on prend la colonne avec le nombre de reads et on la transforme en ligne
CountData = t(file[,6])

#pour ColData on prend deux colonnes : avec le nom et avec le plus ou moins 
ColData = file[,c(1,5)]
#R ne connais pas + et - faut les remplacer par 'plus' et 'minus 
ColData[ColData == '-'] <- 'minus'
ColData[ColData == '+'] <- 'plus'



dds <- DESeqDataSetFromMatrix(countData = CountData,
                              colData = ColData,
                              design= ~Strand)



##Normalement c'est ça qu'il faut faire mais R me met une erreur et me dit de faire le truc d'après à la place 
dds <- DESeq(dds)
resultsNames(dds)

#Ce que me dit de faire R mais le problème c'est qu'après je ne peux pas utiliser les résultats donnés
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)
dispersions(dds) <- mcols(dds)$dispGeneEst



#ne fonctionne pas du coup car pour faire ça il faut faire un DESeq mais R ne veut pas... 
res05 <- results(dds, alpha=0.05)









