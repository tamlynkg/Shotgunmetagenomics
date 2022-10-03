setwd("F:/BMC_Edits/RDA_edit")
library(igraph)
library(ggplot2)
library(reshape2)
rm(list= ls())
library(MASS)
library(vegan)
#read OTU table and metadata
trflp_rda1 <-read.table("otu_table_95.csv",header=T, row.names=1, sep = ",")
trflp_rda1 <-t(trflp_rda1)
env_rda <-read.table("metadata_edit95.csv",header=T, row.names=1, sep = ",")
trflp_rda<-trflp_rda1[row.names(env_rda),]
row.names(trflp_rda)==row.names(env_rda)
#transform otu table
otu.table.t <- decostand(trflp_rda, "hellinger")
assaychem_rda1 <-read.table("meta_95.csv",header=T, row.names=1, sep = ",")
rda_1 <- rda(trflp_rda~., assaychem_rda1)
rda_1
envfit_rda <- envfit(rda_1, assaychem_rda1, perm=999)
envfit_rda

plot(rda_1, type="n")
plot(envfit_rda, p.max=0.08, col= c("black", "red", "blue", "darkgreen"), cex=1)
envfit_rda
env_rda
location <- env_rda[,8]
location
points(rda_1, display = "sites", cex=1, select = which(location=="RF"), pch = 17, col="#DA5724", bg="#DA5724")
points(rda_1, display = "sites", cex=1, select = which(location=="RM"), pch = 19, col="#DA5724", bg="#DA5724")
points(rda_1, display = "sites", cex=1, select = which(location=="UF"), pch = 18, col="#508578", bg="#508578")
points(rda_1, display = "sites", cex=1, select = which(location=="UM"), pch = 15, col="#508578", bg="#508578")

#plot(rda_1)
#add legend
legend('topleft', legend = c("Rural Female", "Rural Male", "Urban Female", "Urban Male"),pch = c(17,19, 18, 15),col = c("#DA5724","#DA5724", "#508578", "#508578"), cex = 1)

