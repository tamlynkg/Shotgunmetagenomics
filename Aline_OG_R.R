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

#function veganotu
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

#Pairwise comparison
library(devtools) #rtools must be installed as well
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis") # to install the plugin pairwise
library(pairwiseAdonis)

### Import Data ####
#-------------------#

setwd("N:/bowiss/RPR/Aline/MetaComp/Data Analyses") # to adapt

#...Metadata
meta <- read.delim(file.choose(), header = TRUE,row.names=1)
str(meta)
rownames(meta)
colnames(meta)
metas=sample_data(meta)

#...bacteria
bact <- read.delim(file.choose(), header = TRUE, row.names=1)
str(bact)
tbact <- as.data.frame(t(bact))
rownames(tbact)
colnames(tbact)
OTUbact=otu_table(tbact,taxa_are_rows=F)

taxobact <- read.delim(file.choose(), header = TRUE, row.names=1)
str(taxobact)
rownames(taxobact)
colnames(taxobact)
taxabact=tax_table(as.matrix(taxobact))

merged_bact = phyloseq(OTUbact, taxabact, metas) #sample names/OTU id must similar when merging data.frames
merged_bact

#...fungi
fung <- read.delim(file.choose(), header = TRUE, row.names=1)
str(fung)
tfung <- as.data.frame(t(fung))
rownames(tfung)
colnames(tfung)
OTUfung=otu_table(tfung,taxa_are_rows=F)

taxofung <- read.delim(file.choose(), header = TRUE, row.names=1)
str(taxofung)
rownames(taxofung)
colnames(taxofung)
taxafung=tax_table(as.matrix(taxofung))

merged_fung = phyloseq(OTUfung, taxafung, meta) #sample names/OTU id must similar when merging data.frames
merged_fung

#...Mesofauna
meso <- read.delim(file.choose(), header = TRUE, row.names=1)
str(meso)
rownames(meso)
colnames(meso)
OTUmeso=otu_table(meso,taxa_are_rows=F)

taxomeso <- read.delim(file.choose(), header = TRUE, row.names=1)
str(taxomeso)
rownames(taxomeso)
colnames(taxomeso)
taxameso=tax_table(as.matrix(taxomeso))

metaMF <- read.delim(file.choose(), header = TRUE,row.names=1)
str(metaMF)
rownames(metaMF)
colnames(metaMF)
metameso=sample_data(metaMF)

merged_meso = phyloseq(OTUmeso, taxameso, metameso) #sample names/OTU id must similar when merging data.frames
merged_meso


### Alphadiversity ####
#---------------------#

#subset samples example
bact1 = subset_samples(merged_bact, site != "Burgdorf") #without samples from site Burgdorf
sample_data(bact1)

#Rarefy to even sequencing depth and calculate indices:
data_to_rarefy=merged_bact   #phyloseq data we need to rarefy.
data_stats=meta
raref_data=rarefy_even_depth(data_to_rarefy,rngseed = T)
str(raref_data)

#Calculate Alpha-diversity Indices
Alphadiv=estimate_richness(raref_data)
Alphadiv=cbind(data_stats,Alphadiv) #combine design data with indices
Alphadiv
write.csv(Alphadiv,"C:/Users/frossard/Desktop/An demo/Bact_Alphadiv.csv",row.names=F) #specify path and file name to save as csv.
str(Alphadiv)
Alphadiv
colnames(Alphadiv)

## summarize alphadiv data for table with diversity indices
colnames(Alphadiv)
Alpha=Alphadiv[,c(1:2,4:5,9:12)] #select indices (columns) to analyze 
Alpha
Alpha.df <- data.frame(Alpha)
sum_alpha <- aggregate(Alphadiv,list(meta$site, meta$Hg), function(x) c(M = mean(x), SE = sd(x)/sqrt(length(x))))
sum_alpha
write.csv(sum_alpha,"C:/Users/frossard/Desktop/An demo/Bact_sumAlpha.csv",row.names=F)


# Scatter plot graph with regression
p=plot_richness(raref_data)
p

p1=plot_richness(raref_data, measures="Observed") #ACE, Shannon, Simpson InvSimpson, Fisher, ..
p1

#barplot example
graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=12),
                    axis.title.y = element_text(size=15, colour="black", angle=90),
                    axis.title.x = element_text(size=15, colour="black", angle=0),
                    plot.title = element_text(size=25))

sumdata <- summarySE(Alphadiv, measurevar="Observed", groupvars=c("siteN","site","Hg"))
sumdata
ggplot(data=sumdata, aes(x=siteN, y=Observed, fill=Hg)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se), width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("Lausanne", "Burgdorf", "Gerlafingen", "Piotta", "Schänis", "Laufen", "Sihlwald"))+
  ggtitle("Bacteria Richness")+
  scale_fill_manual(values = c("white", "grey60", "grey35", "black"))+
  ylab("Richness")+
  graph.theme 

#boxplot example
ggplot(Alphadiv, aes(x = site, y = Observed, fill = Hg)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot() +
  scale_fill_manual(values = c("white", "grey60", "grey35", "black"),
                    limits = c("Hg0", "Hg0.32", "Hg3.2", "Hg32"),
                    label=c("0 ppm", "0.32 ppm", "3.2 ppm", "32 pmm"))+
  scale_x_discrete(label=c("Lausanne", "Burgdorf", "Gerlafingen", "Piotta", "Schänis", "Laufen", "Sihlwald"))+
  #scale_y_log10() + # depending of t he variable take it out
  ggtitle("Bacterial Richness") + 
  ylab("Index")+
  xlab("Sites")+
  #guides(fill=FALSE)+
  graph.theme


# ANOVA for each diversity index: Richness, Shannon, ...
lm1 = lm(Observed ~ site* Hg, data = Alphadiv)
summary.aov(lm1)
TukeyHSD(aov(lm1))
par(mfrow=c(2,3))
plot(lm1)
plot(lm1$fitted,lm1$resid); abline(h=0,lty=2)
hist(lm1$resid)


### Betadiversity ####
#--------------------#

# Ordination with Phyloseq package
#''''''''''''''''''''''''''''''''''#

# standardize data
scaled.data=transform_sample_counts(physeq, function(x) {sqrt(x/sum(x))}) 

# Ordination graph
graph.theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     panel.background = element_rect(colour="black", fill="white"), 
                     axis.line = element_line(colour = "black"),
                     axis.text = element_text(colour="black", size=13),
                     axis.title.y = element_text(size=14, colour="black", angle=90),
                     axis.title.x = element_text(size=14, colour="black"),
                     plot.title = element_text(size=25))

data=scaled.data
ord_type1=ordinate(data,method="PCoA",distance="bray") #change method to obtain different ordinations
ord_type2=ordinate(data,method="NMDS",distance="bray") #change method to obtain different ordinations
ord_type3=ordinate(data,method="RDA",distance="bray") #change method to obtain different ordinations
ord_type4=ordinate(data,method="CAP",distance="bray", formula=otu_table ~ sample_data(data)$site + sample_data(data)$Hg, na.action=na.exclude)
ord_type4b=ordinate(data,method="CAP",distance="bray", formula=otu_table ~ sample_data(data)$site * sample_data(data)$Hg, na.action=na.exclude)

#ordination graph worm part and treat
p1 = plot_ordination(data, ord_type2, color ="Location", shape = "SampleType") +
  scale_colour_manual(values = c("#DA5724", "#508578"),
                      limits = c("Rural", "Urban"),
                      labels = c("Rural", "Urban")) +
  scale_shape_manual(values = c(17, 19, 18, 15),
                     limits = c("RF","RM", "UF", "UM"),
                     label=c("Rural Female","Rural Male", "Urban Female", "Urban Male"))+
  geom_point(size=3) +
  graph.theme
p1

p2 = plot_ordination(data, ord_type3, color ="Location", shape = "SampleType") +
  scale_colour_manual(values = c("#DA5724", "#508578"),
                      limits = c("Rural", "Urban"),
                      labels = c("Rural", "Urban")) +
  scale_shape_manual(values = c(17, 19, 18, 15),
                     limits = c("RF","RM", "UF", "UM"),
                     label=c("Rural Female","Rural Male", "Urban Female", "Urban Male"))+
  geom_point(size=3) +
  graph.theme
p2

ord_type3

#permanova
library("ggplot2")
library("scales")
library("grid")
library(remotes)
library(vegan)
data.veg= veganotu(data) # adapt to vegan package
perma=adonis(data.veg ~ metadata$Location * metadata$SampleType, permutations=999, method="bray") #factors according to design. Use strata=factorXYZ to specify where to permute
perma

pairwise.adonis(data.veg, metadata$SampleType)
pairwise.adonis(data.veg, metadata$Location)
pairwise.adonis(data.veg, metadata$Location : metadata$Gender)


# Ordination using vegan package
#'''''''''''''''''''''''''''''''#

relab = function(x) {sqrt(x/sum(x))} #transform data into relative abundance matrix
scaled.com <- relab(bact)

comx <- bact
str(comx)
envx <- meta
str(envx)

#NMDS: Non metric multidimensional scaling
#' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 
com.nmds <- metaMDS(comx,k=2, trace=1,autotransform=TRUE)
com.nmds


#PCoA: Principal Coordonates Analysis
#' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '
compcoa <- vegdist(comx, "bray")
com.pcoa <- cmdscale(compcoa, k=(nrow(comx)-1), eig=TRUE)
#com.pcoa$values

#RDA: Redundancy analysis
#' ' ' ' ' ' ' ' ' ' ' ' '
com.rda <- rda(comx, trace=1,autotransform=TRUE)
str(com.rda)

#CCA: Constrained correspondance analysis
#' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' '
com.cca1 <- cca(comx ~ moist + temp + day + sub, envx, dist="bray")
com.cca2 <- cca(comx ~ moist * temp * day * sub, envx, dist="bray")
com.cca3 <- cca(comx ~ moist * temp * day + sub, envx, dist="bray")
com.cca3

#CAP: Constrained analysis of proximities => doesnt work right now....
#' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' ' 
com.cap <- capscale(comx ~ moist + temp + day + sub, envx, dist="bray")
summary(com.cap)
com.cap <- capscale(comx ~ moist * temp * day + sub, envx, dist="bray")

#Ordination Plot
#***************   display="sites", type="t")

mycolors <- c("greenyellow", "green4","coral4","orange", "red", "cyan2", "blue4") 
myedge <- c("greenyellow", "green4","coral4","orange", "red", "cyan2", "blue4") 
myshape <- c(21, 24, 22, 23) #define shapes
ordiplot <- plot(com.nmds, display="sites", type="n", main="NMDS Bacteria", choices = c(1, 2)) 
points(ordiplot,"sites", bg=mycolors[envx$site], pch=myshape[envx$Hg], col=myedge[envx$subtemp], cex=1.6)
with(envx, ordiellipse(com.nmds, site, choices=c(1,2), kind = "ehull", conf=0.95, show.groups = "Schänis", col = "red", lty = 5, lwd=2))

#Fitting of environmental variables
#***********************************
ef <- envfit(com.nmds, envx, permu=99999, na.rm=TRUE) #env factors fit
ef

#Permanova on the ordination : adonis
#*************************************
com.ado <- adonis(comx ~  site * Hg, data = envx, method="bray", strata=NULL, perm=9999) # adonis Permanova
com.ado
pa <- pairwise.adonis(comx, meta$site:meta$Hg)# posthoc tests
pa
pa1 <- pairwise.adonis(comx, meta$site)
pa1
pa2 <- pairwise.adonis(comx, meta$Hg)
pa2

### Taxonomy  ####
#----------------#

#Preparation of the taxanomy file
#'''''''''''''''''''''''''''''''''

#upload file
taxo <- read.delim(file.choose(), header = TRUE) # sumary file from pivot table
rownames(taxo)=taxo$sample
str(taxo)
rownames(taxo)
colnames(taxo)

meta <- read.delim(file.choose(), header = TRUE)
rownames(meta)=meta$sample
str(meta)
meta

Tax = cbind(meta, taxo)
str(Tax)
rownames(Tax)
head(Tax)

# mean of the 3 replicates
tax.mean <- aggregate(Tax , by=list(Tax$site, Tax$Hg), FUN = mean)
tax.mean
names(tax.mean)[1] <- "site"
names(tax.mean)[2] <- "Hg"
str(tax.mean)
tax = tax.mean[,-c(3:9)] #take out the unwilling columns
str(tax)
tax

# we need to reshape the matrix into a long matrix
library(reshape2)
tax.long <- melt(tax, id.vars=c("site","Hg"), variable.name = "functional group", value.name = "abundance")
head(tax.long)
names(tax.long)[4] <- "abundance"
names(tax.long)[3] <- "taxalevel"
str(tax.long)
taxlong <- tax.long[order(tax.long$taxalevel),]
head(taxlong)

#Taxonomy plot
#''''''''''''''
library (scales)
library(ggplot2)

taxlong$site <- factor(taxlong$site, levels=c("Lausanne","Burgdorf", "Gerlafingen","Piotta", "Schänis", "Laufen", "Sihlwald"))
taxlong$biome <- factor(taxlong$Hg, levels=c("Hg0","Hg0.32", "Hg3.2", "Hg32"))
str(taxlong)
head(taxlong)

#Example of stacked colum plot Fung Class 1%
ggplot(taxlong, aes(x=Hg, y=abundance, fill=taxalevel)) +
  geom_bar(stat="identity", position = "fill", colour="white", show.legend=TRUE) +
  scale_fill_manual(values=c("indianred1", "indianred3", "indianred4", 
                             "lightblue1", "lightblue2", "lightskyblue", "skyblue3", "royalblue1", "royalblue3", "royalblue4",
                             "peru",
                             "chartreuse3", "chartreuse4",
                             "darkorchid4", "grey50"))+
  scale_y_continuous(labels = percent_format(), name="% phyla")+
  facet_grid(. ~ site )+
  ggtitle("Fung Class 1%cut")

### Indicator Species Analysis ####
#---------------------------------#

#Preparation file
#''''''''''''''''
setwd("D:/WSL/MERCURY/MercuryValais/Data Analyses") # Input files taxo with OTUs integrated (needed for network analyses)

IS_bact <- read.delim(file.choose(), header = TRUE, row.names=1)
str(IS_bact)
rownames(IS_bact)
colnames(IS_bact)
taxab_IS=tax_table(as.matrix(IS_bact))

merged_bact_IS = phyloseq(OTUbact, taxab_IS, metas) #sample names/OTU id must similar when merging data.frames
merged_bact_IS

#List preparation for the results of the ind species analyses
grp = rep(levels(meta$site),each=4) # repeat 4 times because of 4 replicates on the same level = vector cluster
grp

IndSP=list()
PVAL=list()
QVAL=list()
COEFF=list()
COEFF.order=list()
TAX=list()
DF=list()
SB01=list()
SB005=list()
SB001=list()

#Loop for the indicator species analyses
#'''''''''''''''''''''''''''''''''''''''

for (i in 1:4){ #1:4 because 4 treatments (=4 loops)
  
  subset_data = subset_samples(merged_bact_IS,Hg==levels(meta$Hg)[i])  # take subset data, according to the experimental treatment
  
  rel_abund = transform_sample_counts(subset_data, function(x) {sqrt(x/sum(x))}) # transform data with SQRT
  pruned_data = prune_taxa(taxa_sums(subset_data) > 2, subset_data) # pruned data: removing unwanted OTUs such as singleton and doubleton in this case (more than 2)
  Otu_toextract = colnames(veganotu(pruned_data)) #transform phyloseq in vegan matrix
  
  rel_abund_vegan = veganotu(rel_abund)
  rel_abund_ns = rel_abund_vegan[,colnames(rel_abund_vegan) %in% Otu_toextract]
  data_ind <- as.data.frame(rel_abund_ns)
  
  TAX[[i]] = tax_table(pruned_data)
  
  IndSP[[i]]=multipatt(data_ind, cluster= grp, func = "r.g", control = how(nperm=99999))
  
  
  PVAL[[i]]=IndSP[[i]]$sign$p.value
  QVAL[[i]]=p.adjust(PVAL[[i]],method="BH")
  COEFF[[i]]=IndSP[[i]]$str
  COEFF.order[[i]]=COEFF[[i]][order(rownames(COEFF[[i]])),]
  
  DF[[i]]=data.frame(TAX[[i]], COEFF.order[[i]],p.value=PVAL[[i]],q.value=QVAL[[i]])#add Tax[[i]] befor COEFF.order
  
  SB01[[i]]<-subset(DF[[i]], DF[[i]]$q.value<=0.1)
  SB005[[i]]<-subset(DF[[i]], DF[[i]]$q.value<=0.05)
  SB001[[i]]<-subset(DF[[i]], DF[[i]]$q.value<=0.01)
  
}

# loop for saving files

for (i in 1:3){
  
  write.csv(as.data.frame(DF[[i]]),row.names=T,
            paste("D:/WSL/MERCURY/MercuryValais/Data Analyses/IndicSP/",
                  levels(meta$exp.treat)[i],".BACT.allValues.csv",sep=""))
  
  write.csv(IndSP[[i]]$sign,
            paste("D:/WSL/MERCURY/MercuryValais/Data Analyses/IndicSP/",
                  levels(meta$exp.treat)[i],".BACT.sign.csv",sep=""))
  
}

### NETWORK PREPARATION  ####
#----------------------------#

# Co-occurence network
#- - - - - - - - - - - 
#designed function to make a correlation list from correlation matrix
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    sp1 = rownames(cormat)[row(cormat)[ut]],
    sp2 = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# for each treatment and calculate network stats to see the aggregateness of the community.
setwd("D:/WSL/MERCURY/MercuryValais/Data analyses/networks")

str(bact)
rownames(bact)
str(fung)
meta

com <- bact

T1 <- com[c(1:4),]
T2 <- com[c(5:8),]
T3 <- com[c(9:12),]
T4 <- com[c(13:16),]
T1Hg <- com[c(17:20),]
T2Hg <- com[c(21:24),]
T3Hg <- com[c(25:28),]
T4Hg <- com[c(29:32),]
T1c <- com[c(33:36),]
T2c <- com[c(37:40),]
T3c <- com[c(41:44),]
T4c <- com[c(45:48),]

comX <- T4

nr <- rep.int(0,length(comX))
str(nr)
comX.n <- rbind(comX, nr)
rownames(comX.n)
corcom <- comX.n[, colSums(comX.n) >= 10] # keep only OTU with more than 10 sequences
dim(corcom)

# spearman correlation
str(corcom)
library(Hmisc)
res<-rcorr(as.matrix(corcom), type="spearman")
cormat <- flattenCorrMatrix(res$r, res$P)
head(cormat)
str(cormat)

#subsetting the corelation data
cormatsign <- subset(cormat, p <= 0.01)
str(cormatsign)
cordata <- subset(cormatsign, cor <= -0.6 | cor >= 0.6)
head(cordata)
str(cordata)

write.table(cordata, "D:/WSL/MERCURY/MercuryValais/Data analyses/networks/cytoscape_co-occurence_Bact/bact_T4.txt", sep="\t", row.names = F)


# To check the checkerboard combinations (C-score) in a matrix to determine the randomness of the distribution

#C-score: comparison again a null model
library(vegan)
library(bipartite)
# to compare the matrix to a null.model
null.model <- oecosimu(com, bipartite::C.score, "swap", burnin=100, thin=10, statistic="evals", nsimul=100) #where species is you species by sites matrix
print(null.model)
#to get the c-score of the matrix
C.score(com, normalise=T, FUN=mean, na.rm=T)->cscore.speciesN # c-score normaliz between 0 and 1
C.score(com, normalise=F, FUN=mean, na.rm=T)->cscore.speciesS # normal value of the c-score
cscore.speciesN
cscore.speciesS


# Bipartite Network: data = Indicator species data
#- - - - - - - - - -

setwd("D:/WSL/MERCURY/MercuryValais/Data analyses/networks")

com <- read.delim(file.choose(), header = TRUE, row.names=1) # for bipartite network, output of indsp, with corelation to treatments
str(com)
rownames(com)
colnames(com)
#comc = com[,-c(62:71)] #? take out unuseful columns
#str(comc)

com0.3 <- subset(com, T1 > 0.3 | T2 > 0.3 | T3 > 0.3 | T4 > 0.3)
str(com0.3)

com0.3 <- subset(com, T1 < -0.3 | T1 > 0.3 | T2 < -0.3 | T2 > 0.3 | T3 < -0.3 | T3 > 0.3 | T4 < -0.3 | T4 > 0.3)
str(com0.3)

com0.4 <- subset(comc, T1 < -0.4 | T1 > 0.4 | T2 < -0.4 | T2 > 0.4 | T3 < -0.4 | T3 > 0.4 | T4 < -0.4 | T4 > 0.4)
str(com0.4)

com0.5 <- subset(comc, T1 < -0.5 | T1 > 0.5 | T2 < -0.5 | T2 > 0.5 | T3 < -0.5 | T3 > 0.5 | T4 < -0.5 | T4 > 0.5)
str(com0.5)

comX <- com0.3


library(data.table)
colnames(comX)
meltcom = comX[,c(9:12)] #select column for edge table
head(meltcom)
setDT(meltcom, keep.rownames = TRUE)[]# to make the first column as value
colnames(meltcom)[1] <- "OTU"
head(meltcom)

longnet <- melt(meltcom)
head(longnet)
tail(longnet)
str(longnet)
#long <- subset(longnet, value < -0.5 | value > 0.5) # !! change corr. factor!! subsetting to have only links (edge)> 0.5
long <- subset(longnet, value > 0.3) 
str(long)
head(long)
long[order(long$OTU), ]
head(long)

#  !!!Change name!!!
write.csv(long,"D:/WSL/MERCURY/MercuryValais/Data analyses/networks/Bact_field_p0.05_c0.3.csv",row.names=F)

#Node table
colnames(comX)
node <- comX[,c(2:61)]# node table
setDT(node, keep.rownames = TRUE)[]# to make the first column as value
colnames(node)[1] <- "OTU"
str(node)
head(node)
write.csv(node,"D:/WSL/MERCURY/MercuryValais/Data analyses/networks/Bact_field_p0.05_c0.3_NODE.csv",row.names=F)


# Taxonomy network
#- - - - - - - - #

# See special script

### Network Statistics ####
#--------------------------#

setwd("D:/WSL/MERCURY/MercuryValais/Data Analyses/networks/cytoscape_co-occurence_Fung")

# Putting all the data together
net <- read.csv(file.choose(), header = TRUE, row.names=1)
str(net)
dim(net)
netname <- rep("T2Hg",(nrow(net)))  # Add a column to identify the network: change name and number of repetition
str(netname)
type <- rep("Hg", nrow(net)) # type of sample
str(type)
Hglevel <- rep("moderate",nrow(net)) #Hg level
str(Hglevel)
netc <- cbind(netname, type, Hglevel, net)
head(netc)
str(netc)
netT2Hg <- netc #load all networks in the R env. names = T1, T2, T3, T4, T11Hg, T2Hg, T3Hg, T4Hg, T1c, T2c, T3c, T4c
write.table(netT2Hg, "D:/WSL/MERCURY/MercuryValais/Data analyses/networks/cytoscape_co-occurence_Fung/NetStatT2Hg.txt", sep="\t", row.names = F)

head(netT4Hg)

NetStatHgValais <- rbind(netT1, netT2, netT3, netT4, netT1Hg, netT2Hg, netT3Hg, netT4Hg, netT1c, netT2c, netT3c, netT4c)
head(NetStatHgValais)
str(NetStatHgValais)
colnames(NetStatHgValais)
#names(NetStatHgValais)[1] <- "netname"

write.table(NetStatHgValais, "D:/WSL/MERCURY/MercuryValais/Data analyses/networks/cytoscape_co-occurence_Fung/StatNetHgValais.txt", sep="\t", row.names = F)

# summary of the data
library(dplyr)
library(magrittr)
library(ggplot2)

stderr <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

colnames(NetStatHgValais)

sumstat <- NetStatHgValais %>%
  group_by(netname) %>%
  summarise(nbConn.mean   = mean(NumberOfDirectedEdges),
            nbConn.stderr = stderr(NumberOfDirectedEdges),
            Topol.mean   = mean(TopologicalCoefficient),
            Topol.stderr = stderr(TopologicalCoefficient),
            AvshortPath.mean   = mean(AverageShortestPathLength),
            AvshortPath.stderr = stderr(AverageShortestPathLength),
            Clust.mean   = mean(ClusteringCoefficient),
            Clust.stderr = stderr(ClusteringCoefficient),
            ClosCentral.mean   = mean(ClosenessCentrality),
            ClosCentral.stderr = stderr(ClosenessCentrality),
            Eccentr.mean   = mean(Eccentricity),
            Eccentr.stderr = stderr(Eccentricity)
  )
str(sumstat)        
sumstat
write.table(sumstat, "D:/WSL/MERCURY/MercuryValais/Data analyses/networks/cytoscape_co-occurence_Fung/sumstatCo-occurNet.txt", sep="\t", row.names = F)

# graphs
library(plyr)
library(Rmisc)

setwd("D:/WSL/MERCURY/MercuryValais/Data Analyses/networks/cytoscape_co-occurence_Fung")

NetStatHgValais <- read.delim(file.choose(), header = TRUE)
str(NetStatHgValais)

graph.theme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                    panel.background = element_rect(colour="black", fill="white"), 
                    axis.line = element_line(colour = "black"),
                    axis.text = element_text(colour="black", size=15),
                    axis.title.y = element_text(size=20, colour="black", angle=90),
                    axis.title.x = element_blank(),
                    plot.title = element_text(size=25))

#variables:
#NumberOfDirectedEdges
#TopologicalCoefficient
#AverageShortestPathLength
#ClusteringCoefficient
#ClosenessCentrality
#Eccentricity

str(NetStatHgValais)
colnames(NetStatHgValais)

###
sum.data <- summarySE(NetStatHgValais, measurevar="NumberOfDirectedEdges", groupvars=c("Hglevel","type"))
str(sum.data)

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=NumberOfDirectedEdges, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=NumberOfDirectedEdges-se, ymax=NumberOfDirectedEdges+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(6,14))+
  ylab("Number of Edges")+
  ggtitle("Edges")+
  graph.theme 

###
sum.data <- summarySE(NetStatHgValais, measurevar="TopologicalCoefficient", groupvars=c("Hglevel","type"))
sum.data

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=TopologicalCoefficient, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=TopologicalCoefficient-se, ymax=TopologicalCoefficient+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(0.6,0.9))+
  ylab("Coefficient")+
  ggtitle("Topological Coefficient")+
  graph.theme  

###
sum.data <- summarySE(NetStatHgValais, measurevar="AverageShortestPathLength", groupvars=c("Hglevel","type"))
sum.data

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=AverageShortestPathLength, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=AverageShortestPathLength-se, ymax=AverageShortestPathLength+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(1,5))+
  ylab("number of connection")+
  ggtitle("Shortest Pathlength")+
  graph.theme  

###
sum.data <- summarySE(NetStatHgValais, measurevar="ClusteringCoefficient", groupvars=c("Hglevel","type"))
sum.data

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=ClusteringCoefficient, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=ClusteringCoefficient-se, ymax=ClusteringCoefficient+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(0.8,0.95))+
  ylab("Coefficient")+
  ggtitle("Clustering Coefficient")+
  graph.theme 

###
sum.data <- summarySE(NetStatHgValais, measurevar="ClosenessCentrality", groupvars=c("Hglevel","type"))
sum.data

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=ClosenessCentrality, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=ClosenessCentrality-se, ymax=ClosenessCentrality+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(0.4,0.8))+
  ylab("Coefficient")+
  ggtitle("Centrality")+
  graph.theme 

###
sum.data <- summarySE(NetStatHgValais, measurevar="Eccentricity", groupvars=c("Hglevel","type"))
sum.data

sum.data$type <- factor(sum.data$type, levels = c("field", "Hg", "ctrl"))
sum.data$type

ggplot(data=sum.data, aes(x=Hglevel, y=Eccentricity, fill=type)) + 
  geom_bar(position=position_dodge(), stat = "identity", colour="black") +
  geom_errorbar(aes(ymin=Eccentricity-se, ymax=Eccentricity+se), 
                width=0.25, position=position_dodge(0.9))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_fill_manual(values = c("black", "grey60", "white"), labels=c("Field", "Hg", "Control"))+
  coord_cartesian(ylim=c(2,10))+
  ylab("Coefficient")+
  ggtitle("Eccentricity")+
  graph.theme 


# boxplot

ggplot(sumdata, aes(factor(Hg), AvshortPath)) + 
  geom_boxplot(aes(fill=factor(Domain)), colour="black", size=0.5, outlier.size=1)+ 
  scale_fill_manual(values = c("grey60", "black"))+
  ggtitle("Averag shortest Pathway") + 
  ylab("Shortest pahtway")+
  guides(fill=FALSE)+
  graph.theme

# Boxplot for the entire database

str(NetStatHgValais)
colnames(NetStatHgValais)


ggplot(NetStatHgValais, aes(x = Hglevel, y = NumberOfDirectedEdges, fill = type)) + 
  stat_boxplot(geom ='errorbar')+
  geom_boxplot() +
  scale_fill_manual(values = c("black", "grey60", "white"),
                    limits = c("field", "ctrl", "Hg"),
                    label=c("Field", "control", "Hg"))+
  scale_x_discrete(label=c("high", "moderate", "low", "natural"))+
  scale_y_continuous(limits = c(20, 65)) + 
  ggtitle("Edges") + 
  ylab("Number of Edges")+
  xlab("Level of Hg contaminated sites")+
  guides(fill=FALSE)+
  graph.theme


# ANOVA

setwd("D:/WSL/MERCURY/MercuryValais/Data Analyses/networks")

statnet <- read.delim(file.choose(), header = TRUE) # dataset from stat network cytoscape
str(statnet)

#select variables:NumberOfDirectedEdges,TopologicalCoefficient,AverageShortestPathLength,ClusteringCoefficient,ClosenessCentrality,Eccentricity)
stat <- subset(statnet, select=c(NumberOfDirectedEdges,TopologicalCoefficient,
                                 AverageShortestPathLength,ClusteringCoefficient,
                                 ClosenessCentrality,Eccentricity))
str(stat)

lgstatnet <- log(stat+1)
str(lgstatnet)

var <- subset(statnet, select=c(type, Hglevel))

trstatnet <- cbind(var,lgstatnet) 
str(trstatnet)

x <- trstatnet$Eccentricity
str(x)

lm1 = lm(x ~ Hglevel * type, data = trstatnet)
summary.aov(lm1)
TukeyHSD(aov(lm1))
par(mfrow=c(2,2))
plot(lm1)
qqnorm(lm1$resid); qqline(lm1$resid)
plot(lm1$fitted,lm1$resid); abline(h=0,lty=2)
hist(lm1$resid)



