### Diversity Dataset Analyses ###
##################################

### Load libraries ####
#---------------------#
library(vegan)
library(labdsv)
library(indicspecies) 
library(ggplot2)
library(multcomp)
library(multcompView)
library(reshape)
library(phyloseq)
library(agricolae)
library(RVAideMemoire)
library(gplots)
library (MASS)
library (rgl)
library (lattice)
library(ape)
library(plyr)
library(Rmisc)
library(tidyverse)
library(ggpubr)
library(pairwiseAdonis)

#function veganotu
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}


### Import Data for 95 samples####
#-------------------#
setwd("F:/")#CHONGO

otu_table = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#otu_table_OG_95
taxonomy = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#Taxa
metadata = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#metadata_edit95
otu_table = as.matrix(otu_table)
taxonomy = as.matrix(taxonomy)
class(otu_table)
class(taxonomy)
#Recognize your files as phyloseq objects 
otu = otu_table(otu_table, taxa_are_rows = TRUE)
tax = tax_table(taxonomy)
sampledata = sample_data(metadata)
#Merging three tables into one phyloseq
physeq = phyloseq(otu, tax, sampledata)
#check the object physeq
physeq

### Alphadiversity ####
#---------------------#

#Rarefy to even sequencing depth and calculate indices:
data_to_rarefy=physeq   #phyloseq data we need to rarefy.
data_stats=sampledata
raref_data=rarefy_even_depth(data_to_rarefy,rngseed = T)
str(raref_data)

#Calculate Alpha-diversity Indices
Alphadiv=estimate_richness(physeq)
data_stats=sampledata
Alphadiv95=cbind(data_stats,Alphadiv) #combine design data with indices
Alphadiv95
write.csv(Alphadiv95,"F:/Alphadiv95.csv",row.names=F) #specify path and file name to save as csv.
str(Alphadiv95)
Alphadiv95
colnames(Alphadiv95)

## summarize alphadiv data for table with diversity indices
Alpha=Alphadiv95[,c(1:2,5:22)] #select indices (columns) to analyze 
Alpha
Alpha.df <- data.frame(Alpha)
sum_alpha1 <- aggregate(Alphadiv1,list(meta$Location, meta$Gender, meta$SampleType, meta$Age_group), function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
sum_alpha1
write.csv(sum_alpha1,"F:/Fung_sumAlpha1.csv",row.names=F)

#barplot example
graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=10),
                    axis.title.y = element_text(size=13, colour="black", angle=90),
                    axis.title.x = element_text(size=13, colour="black", angle=0),
plot.title = element_text(size=15))

#Observed grouped by location
Alphadiv95 <- read.csv(file.choose()) 
sumdata1 <- summarySE(Alphadiv95, measurevar="Observed", groupvars=c("Location"))
sumdata1
P1 = ggplot(data=sumdata1, aes(x=Location, y=Observed, fill=Location)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("Rural", "Urban"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#508578","#DA5724"))+
  ylab("Observed species richness")+
  graph.theme
P1
#Wilcoxon test
res1 <- wilcox.test(Observed ~ Location, data = Alphadiv95,
                    exact = FALSE)
res1
#kruskal wallis test
res2<- kruskal.test(Observed ~ Location, data = Alphadiv95)
res2
#Observed by gender
sumdata2 <- summarySE(Alphadiv95, measurevar="Observed", groupvars=c("Gender"))
sumdata2
P2 = ggplot(data=sumdata2, aes(x=Gender, y=Observed, fill=Gender)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("Female", "Male"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578"))+
  ylab("Observed species richness")+
  graph.theme
P2

res3 <- wilcox.test(Observed ~ Gender, data = Alphadiv95, exact = FALSE)
res3
res4<- kruskal.test(Observed ~ Gender, data = Alphadiv95)
res4
#Observed grouped by SampleType
sumdata3 <- summarySE(Alphadiv95, measurevar="Observed", groupvars=c("SampleType"))
sumdata3
P3 = ggplot(data=sumdata3, aes(x=SampleType, y=Observed, fill=SampleType)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("RF", "RM", "UF", "UM"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578", "#5F7FC7", "#CD9BCD"))+
  ylab("Observed species richness")+
  graph.theme
P3

res5<- kruskal.test(Observed ~ SampleType, data = Alphadiv95)
res5
#Observed by Age_group
sumdata4 <- summarySE(Alphadiv95, measurevar="Observed", groupvars=c("Age_group"))
sumdata4
P4 = ggplot(data=sumdata4, aes(x=Age_group, y=Observed, fill=Age_group)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("MAA", "OA", "YA"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578", "#CD9BCD"))+
  ylab("Observed species richness")+
  graph.theme
P4

res6<- kruskal.test(Observed ~ Age_group, data = Alphadiv95)
res6
library(gridExtra)
grid.arrange(P1, P2, P3, P4, ncol=4)

######################################For 100 samples#########

otu_table = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#OTU_table
taxonomy = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#Taxa
metadata = read.table(file.choose(), header = TRUE, row.names = 1, sep = "\t")#metadata
otu_table = as.matrix(otu_table)
taxonomy = as.matrix(taxonomy)
class(otu_table)
class(taxonomy)
#Recognize your files as phyloseq objects 
otu = otu_table(otu_table, taxa_are_rows = TRUE)
tax = tax_table(taxonomy)
sampledata = sample_data(metadata)
#Merging three tables into one phyloseq
physeq = phyloseq(otu, tax, sampledata)
#check the object physeq
physeq

### Alphadiversity ####
#---------------------#
data_stats=sampledata

#Calculate Alpha-diversity Indices
Alphadiv=estimate_richness(physeq)
data_stats=sampledata
Alphadiv100=cbind(data_stats,Alphadiv) #combine design data with indices
Alphadiv100
write.csv(Alphadiv100,"F:/BMC_Edits/Alphadiv100.csv",row.names=F) #specify path and file name to save as csv.
str(Alphadiv100)
Alphadiv100
colnames(Alphadiv100)

#barplot example
graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=10),
                    axis.title.y = element_text(size=13, colour="black", angle=90),
                    axis.title.x = element_text(size=13, colour="black", angle=0),
                    plot.title = element_text(size=15))

#Observed grouped by location
Alphadiv100 <- read.csv(file.choose()) 
sumdata1 <- summarySE(Alphadiv100, measurevar="Observed", groupvars=c("Location"))
sumdata1
P1 = ggplot(data=sumdata1, aes(x=Location, y=Observed, fill=Location)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("Rural", "Urban"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#508578","#DA5724"))+
  ylab("Observed species richness")+
  graph.theme
P1
#Wilcoxon test
res1 <- wilcox.test(Observed ~ Location, data = Alphadiv100,
                    exact = FALSE)
res1
#kruskal wallis test
res2<- kruskal.test(Observed ~ Location, data = Alphadiv100)
res2
#Observed by gender
sumdata2 <- summarySE(Alphadiv100, measurevar="Observed", groupvars=c("Gender"))
sumdata2
P2 = ggplot(data=sumdata2, aes(x=Gender, y=Observed, fill=Gender)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("Female", "Male"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578"))+
  ylab("Observed species richness")+
  graph.theme
P2

res3 <- wilcox.test(Observed ~ Gender, data = Alphadiv100, exact = FALSE)
res3
res4<- kruskal.test(Observed ~ Gender, data = Alphadiv100)
res4
#Observed grouped by SampleType
sumdata3 <- summarySE(Alphadiv100, measurevar="Observed", groupvars=c("SampleType"))
sumdata3
P3 = ggplot(data=sumdata3, aes(x=SampleType, y=Observed, fill=SampleType)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("RF", "RM", "UF", "UM"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578", "#5F7FC7", "#CD9BCD"))+
  ylab("Observed species richness")+
  graph.theme
P3

res5<- kruskal.test(Observed ~ SampleType, data = Alphadiv100)
res5
#Observed by Age_group
sumdata4 <- summarySE(Alphadiv100, measurevar="Observed", groupvars=c("Age_group"))
sumdata4
P4 = ggplot(data=sumdata4, aes(x=Age_group, y=Observed, fill=Age_group)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(1.5))+
  scale_x_discrete(label=c("MAA", "OA", "YA"))+
  ggtitle("Observed")+
  scale_fill_manual(values = c("#DA5724", "#508578", "#CD9BCD"))+
  ylab("Observed species richness")+
  graph.theme
P4

res6<- kruskal.test(Observed ~ Age_group, data = Alphadiv100)
res6
library(gridExtra)
grid.arrange(P1, P2, P3, P4, ncol=4)
