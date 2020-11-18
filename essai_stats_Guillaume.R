
## Installer les bonnes librairies :

install.packages("htmltools")
install.packages("ggplot2")
library(htmltools)
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("questionr")

# Charger les librairies :
library( "DESeq2")
library(ggplot2)

########################################################################################################
## La fonction DESeqDataSetFromMatrix prends en entrée :                                               #
#                                                                                                      #
#     - une matrice de nombres qui sont les counts                                                     #
#     - le nom des colonnes au format data.frame, pour pouvoir faire les différences selon le nom      #
#     - un shéma d'analyse, et c'est là que j'ai pas compris ce qu'on devait mettre                    #
########################################################################################################

###### Récupérer pour chaque sample les counts :

test1 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628582.counts")
test2 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628583.counts")
test3 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628584.counts")
test4 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628585.counts")
test5 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628586.counts")
test6 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628587.counts")
test7 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628588.counts")
test8 = read.table("C:/Users/Guips/Desktop/countfiles/SRR628589.counts")

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
labels2 = data.frame(labels)
rowlabels = test1[-1,1]

###### ajout des labels des colonnes et lignes de la matrice :
rownames(essai)=rowlabels
colnames(essai) = labels

###### regarder à quoi ça ressemble :
head(essai)


###### Etude : ce qui va pas, c'est ce qu'il y a après la petite vague parce que je sais pas ce que ça fait ! ^^
dds = DESeqDataSetFromMatrix(countData = essai, colData = labels2, ~1)
dds2 = DESeq(dds)
resultats = results(dds2)

head(resultats)
