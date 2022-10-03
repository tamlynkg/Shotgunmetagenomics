library(vegan)

# read the table
soil_chem <- read.delim("metadata.txt", header = T, row.names = 1)
# attach the table for easy access

attach(soil_chem)

#check the col names
names(soil_chem)

as.numeric(Inflammation_Status)
##check the data by bnoxplot
boxplot(Inflammation_Status~BV.Status, data = soil_chem, xlab = "BV Status", ylab = "Inflammation Status")

?boxplot
##anova of moisture by carbon
anova_carbon <-aov(Inflammation_Status~BV.Status, data = soil_chem)
# Stats summary
summary(anova_carbon)
#tukey HSD post hoc test
tukey.plot.aov <- TukeyHSD(anova_carbon)

plot(tukey.plot.aov, las = 0.5, cex.axis = 1)
??las
?plot
# AnOVA and Tukey HSD for pH
anova_pH <-aov(pH~Moisture, data = soil_chem)
summary(anova_pH)
tukey_pH <- TukeyHSD(anova_pH)
