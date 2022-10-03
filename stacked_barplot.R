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
#Removing unwanted samples if any
physeq_95 <- subset_samples(physeq,sample_names(physeq) != "T45" & sample_names(physeq) != "T46" & sample_names(physeq) != "T47" & sample_names(physeq) != "T49" & sample_names(physeq) != "T50")
physeq_95
### Alphadiversity ####
#---------------------#


#Calculate Alpha-diversity Indices
Alphadiv3=estimate_richness(physeq_95)
data_stats=sampledata
Alphadiv3=cbind(data_stats,Alphadiv3) #combine design data with indices
Alphadiv3
write.csv(Alphadiv3,"F:/BMC_Edits/Alphadiv3.csv",row.names=F) #specify path and file name to save as csv.
str(Alphadiv3)
Alphadiv3
colnames(Alphadiv3)

#barplot example
graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=10),
                    axis.title.y = element_text(size=13, colour="black", angle=90),
                    axis.title.x = element_text(size=13, colour="black", angle=0),
                    plot.title = element_text(size=15))

#Observed grouped by location
Alphadiv3 <- read.csv(file.choose()) 
sumdata1 <- summarySE(Alphadiv3, measurevar="Observed", groupvars=c("Location"))
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
res1 <- wilcox.test(Observed ~ Location, data = Alphadiv3,
                    exact = FALSE)
res1
#kruskal wallis test
res2<- kruskal.test(Observed ~ Location, data = Alphadiv3)
res2

#Observed by gender
sumdata2 <- summarySE(Alphadiv3, measurevar="Observed", groupvars=c("Gender"))
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

res3 <- wilcox.test(Observed ~ Gender, data = Alphadiv, exact = FALSE)
res3
res4<- kruskal.test(Observed ~ Gender, data = Alphadiv)
res4
#Observed grouped by SampleType
sumdata3 <- summarySE(Alphadiv3, measurevar="Observed", groupvars=c("SampleType"))
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

res5<- kruskal.test(Observed ~ SampleType, data = Alphadiv3)
res5
#Observed by Age_group
sumdata4 <- summarySE(Alphadiv3, measurevar="Observed", groupvars=c("Age_group"))
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

res6<- kruskal.test(Observed ~ Age_group, data = Alphadiv3)
res6
library(gridExtra)
grid.arrange(P1, P3, P4, ncol=3)

#Stacked barplots for phylum
#****************************
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(devtools)
prune_phylum2 <- physeq_95 %>%
  tax_glom(taxrank="Phylum") %>%                            # agglomerate at phylum level
  transform_sample_counts(function(x) {x*100/sum(x)}) %>%   # Transform to rel. abundance
  psmelt() %>%                                              # Melt to long format
  filter(Abundance > 0.1) %>%                                 # Filter out low abundance taxa
  arrange(Phylum)
write.csv(prune_phylum2, "prune_phylum2.csv")
phylum_colors <- c(phylum_colors <- c("#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#999999", "#E69F00", "#56B4E9", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

prune_phylum2$Phylum <- factor(prune_phylum2$Phylum, levels = unique(prune_phylum2$Phylum)) #preserve order of the phyla
graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=10),
                    axis.title.y = element_text(size=13, colour="black", angle=90),
                    axis.title.x = element_text(size=13, colour="black", angle=0),
                    plot.title = element_text(size=15))
p1 = ggplot(prune_phylum2, aes(x=Location, y=Abundance, fill=Phylum)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=phylum_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x") + graph.theme
p1
#Wilcoxon test
resp1 <- wilcox.test(Abundance ~ Location, data = prune_phylum2,
                    exact = FALSE)
resp1


p2 = ggplot(prune_phylum2, aes(x=Age_group, y=Abundance, fill=Phylum)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=phylum_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Phyla > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x")
p2
#kruskal wallis test
resp2<- kruskal.test(Abundance ~ Age_group, data = prune_phylum2)
resp2


#######Class
prune_class2 <- physeq_95 %>%
  tax_glom(taxrank="Class") %>%                            # agglomerate at phylum level
  transform_sample_counts(function(x) x/sum(x)) %>%   # Transform to rel. abundance
  psmelt() %>%                                              # Melt to long format
  filter(Abundance > 0.1) %>%                                 # Filter out low abundance taxa
  arrange(Class)
write.csv(prune_class2, "prune_class.csv")
phylum_colors <- c(phylum_colors <- c("#CD9BCD", "#508578", "#AD6F3B", "#673770", "#D14285", "#DA5724", "#652926", "#C84248", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#999999", "#E69F00", "#56B4E9", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

prune_class2$Class <- factor(prune_class2$Class, levels = unique(prune_class2$Class)) #preserve order of the phyla

p3 = ggplot(prune_class2, aes(x=Location, y=Abundance, fill=Class)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=phylum_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x") + graph.theme 
p3
resp3 <- wilcox.test(Abundance ~ Location, data = prune_class2,
                    exact = FALSE)
resp3
p4 = ggplot(prune_class, aes(x=SampleType, y=Abundance, fill=Class)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=phylum_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Class > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x")
p4
#kruskal wallis test
resp4<- kruskal.test(Abundance ~ Age_group, data = prune_class2)
resp4

#######genus
prune_genus <- physeq_95 %>%
  tax_glom(taxrank="Genus") %>%                            # agglomerate at phylum level
  transform_sample_counts(function(x) x/sum(x)) %>%   # Transform to rel. abundance
  psmelt() %>%                                              # Melt to long format
  filter(Abundance > 0.1) %>%                                 # Filter out low abundance taxa
  arrange(Genus)
write.csv(prune_genus, "prune_genus.csv")
genus_colors <- c(genus_colors <- c("#DA5724", "#508578", "#CD9BCD","#C59B76", "#AD6F3B", "#673770","#C0C781", "#D14285", "#652926", "#C84248",
                                    "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "#999999", "#E69F00", "#56B4E9",
                                    "#999999", "#E69F00", "#56B4E9", "#009E73", "#3277a8", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                                    "#32a852", "#D77363", "#C03221", "#F2D0A4", "#545E75", "#3F826D", "#B279A7", "#ED6A5E", "#FFAF87"))

prune_genus$Genus <- factor(prune_genus$Genus, levels = unique(prune_genus$Genus)) #preserve order of the phyla

p5 = ggplot(prune_genus, aes(x=Location, y=Abundance, fill=Genus)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=genus_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x")
p5
resp5 <- wilcox.test(Abundance ~ Location, data = prune_genus,
                     exact = FALSE)
resp5
p6 = ggplot(prune_genus, aes(x=SampleType, y=Abundance, fill=Genus)) + #x= (change to samples for individual columns)
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=genus_colors) +
  scale_y_continuous(labels=scales::percent)+
  theme(axis.title.x = element_blank()) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genus > 0.1%) \n") +
  theme_bw() + facet_grid(~Location, scales = "free_x")
p6
#kruskal wallis test
resp6<- kruskal.test(Abundance ~ Age_group, data = prune_genus)
resp6
grid.arrange(p1, p3, p5, ncol=3)
grid.arrange(p2, p4, ncol=2)
