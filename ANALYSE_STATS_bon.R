## Installer les bonnes librairies :

#install.packages("htmltools")
#install.packages("ggplot2")
#library(htmltools)
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install.packages("questionr")

# Charger les librairies :
library( "DESeq2")
library(ggplot2)

########################################################################################################
# La fonction DESeqDataSetFromMatrix prends en entrée :                                                #
#                                                                                                      #
#     - une matrice de nombres qui sont les counts                                                     #
#     - le nom des colonnes au format data.frame, pour pouvoir faire les différences selon le nom      #
#     - un shéma d'analyse, et c'est là que j'ai pas compris ce qu'on devait mettre                    #
########################################################################################################

###### Récupérer pour chaque sample les counts :

test1 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628582.counts")
test2 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628583.counts")
test3 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628584.counts")
test4 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628585.counts")
test5 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628586.counts")
test6 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628587.counts")
test7 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628588.counts")
test8 = read.table("C:/Users/Administrator/Desktop/AMI2B/Hackathon/countfiles/SRR628589.counts")

##############################################################################################################################################
# les infos sur les counts sont dans la dernière colonne, donc on les réunit dans une seule matrice après les avoir transformés en numeric : #
##############################################################################################################################################

###### assemblage des dernièeres lignes de chaque fichier en les tranformant en numeric :
essai = data.frame(as.numeric(test1[-1,7]), as.numeric(test2[-1,7]), as.numeric(test3[-1,7])
                   , as.numeric(test4[-1,7]), as.numeric(test5[-1,7]), as.numeric(test6[-1,7]), 
                   as.numeric(test7[-1,7]), as.numeric(test8[-1,7]))

###### passage au format matrix:
essai = as.matrix(essai)

###### récupération au format data.frame des labels :
labels = c("SRR628582", "SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")
mutations = c("M", "M", "WT", "WT", "WT", "WT", "WT", "M")
labels2 = data.frame(labels, mutations)
rowlabels = test1[-1,1]

###### ajout des labels des colonnes et lignes de la matrice :
rownames(essai)=rowlabels
colnames(essai) = labels

###### regarder à quoi ça ressemble :
head(essai)


###### Etude 
dds = DESeqDataSetFromMatrix(countData = essai, colData = labels2, ~mutations)
dds2 = DESeq(dds)
resultats = results(dds2)

head(resultats)

#par(mfrow=c(1,1))
#with(resultats, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-8,8)))
#with(subset(resultats, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
#with(subset(resultats, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

summary(resultats)

##Pour installer le package EnhancedVolcano
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')


##Comme on ne veut pas afficher de noms car cela rend le graphiqu emoins lisible
#On veut lui dire de n'afficher les noms que des gènes étant très à gauche ou très à droite du graphe
keyvals <- ifelse(resultats$log2FoldChange < -10, 'petit',
 ifelse(resultats$log2FoldChange > 10, 'grand','normal'))

#library(EnhancedVolcano)
EnhancedVolcano(resultats,
                lab = rownames(resultats),
                x = 'log2FoldChange', selectLab = rownames(resultats)[which(keyvals %in% c('petit', 'grand'))],
                y = 'padj', 
                title = 'WT VS muté avec p-value ajustée', xlim = c(-10,10), ylim = c(-1,20),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0, drawConnectors = TRUE,widthConnectors = 0.75,
                labSize = 6.0)



########################################
# Analyse en composantes principales : #
########################################


#install.packages("FactoMineR")
#install.packages("factoextra")
library("FactoMineR")
library("factoextra")

essai2 = t(essai)
res_pca2 = PCA(X = essai2, graph = FALSE)
# --> regarder le PCA graph of individuals : on a bien 5 group�s et 3 pas group�s ! Et c'est les bons (mut�s) qui sont pas group�s!
fviz_pca_ind (res_pca2,
              repel = TRUE, # �vite le chevauchement de texte, 
              col.ind = mutations, # colorer by groups
              palette = c("#00AFBB", "#E7B800"),
              legend.title = "Groups",
              mean.point = FALSE
)

