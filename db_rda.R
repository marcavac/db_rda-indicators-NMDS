library(devtools)
library(lattice)
library(multcomp)
library(sandwich)
library(ggplot2)
library(vegan)
library(CCA)
library(tidyr)
library(CCP)
library(GGally)
#function VIF is in library car
library(car)
library(caret)
library(dplyr)
library(phyloseq)
library(data.table)
library(gridExtra)
library(desiderata)
library(FactoMineR)
library(MASS)
library(ellipse)
library(ade4)
library(esquisse)
library(desiderata)
library(Rfast)
#-----------------------------------------RDA-- using data from Liam's permafrost project----------------------------------
#steps:
#read in asv and metadata
#remove NAs/text variables in metadata
#match row names between asv and metadata
# scale env data
#run db rda model
#check collinearity
#run backwards or forwards stepwise model
#run stats
#plot

#----------------------------------------------------------------------------
#methanogens-read in text
meth<-read.csv("methanogen_otu_taxa.csv", row.names=1)
#remove taxa columns:
meth<-meth[,-c(1,48:52)]
tmeth<-t(meth)
#"normalize" dataset
tmeth<- decostand(tmeth,"hellinger")
#read in metadata:
metadata2<-read.csv("all_OTU_NMDS.csv", row.names=1)
##get rid of unwanted columnds
metadata2<-metadata2[,-c(1:3,5:10, 12, 13)]
#make rownames same between meta and asv tables
metadata2<-metadata2[rownames(tmeth),] #make sure data makes sense/ has no NAs.
#remova all but carbon *removal after testing VIF below)
metadata2<-metadata2[,-c(9:16)]
#remove SUVA
metadata2<-metadata2[,-5]
metadata2<-na.omit(metadata2)
#match rownames in both datasets
tmeth<-tmeth[rownames(metadata2),] # make sure data makes sense- may only need to this one or the above, likely not both...
#----------------------------scale env var for select things
#select specific metadata variables
metadata2<-metadata[,c(6,8,14,15,17,18,31,20)]
metadata2<-na.omit(metadata2)
#select categorical variables
#metadata4<-metadata2[,c(3,4,11,10,12,14,15,16,17,28,6)]
#metadata4<-na.omit(metadata4)
#scale variables
env.z <- scale(metadata2)
env.z<-as.data.frame(env.z)
#remove depth if using distance to water table ------ additional processing of my own dataset- ignore ****
#env.z<-env.z[,-4]
#remove water  table if using distance to water table
#env.z<-env.z[,-1]
#remove SUVA
#metadata2<-metadata4[,-8]
#remove SUVA, methanogen:
#env.z<-env.z[,-4]
#remove isotopic signatures
#env.z<-env.z[,-c(6,7)]
# the above breaks b/c we have a categorical factor in env 
# vegan requires that we write out each term if we are not going to 
# convert the factor to a dummy matrix 
set.seed(1234)
rda_tree = rda(tmeth ~ . , data=env.z, scale=T)
summary(rda_tree)
R2<-RsquareAdj(rda_tree)
#R2 adj measures amount of unbiased explained var
(R2adj <- RsquareAdj(rda_tree)$adj.r.squared)# here, model explains 17% of variatio
#plot(rda_tree, type='n', scaling=1)
#orditorp(rda_tree, display='sites', cex=0.5, scaling=1, col='blue')
#text(rda_tree, display='cn', col='red')
rda<-rda(meth)
biplot(rda, choice=c(1,2), display="sites")
#autoplot(rda_tree, arrows = TRUE)
#alternative RDA plotting:, global RDA
#https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
par(mar=c(4,4,2,2))
plot(rda_tree, xlab="RDA1 (15.11 %)", ylab="RDA2 (13.70 %)", 
     display=c("cn", "lc", "sp"), type="n", xlim=c(-1.3,1.3))
sites.sc <- scores(rda_tree, choices=1:2, scaling=2, display="lc")
points(sites.sc, pch=1, cex=5)
text(rda_tree, display = "sites")
#shows species as red arrows
#spe.sc <- scores(rda_tree, choices=1:2, scaling=1, display="sp")
#arrows(0,0, spe.sc[,1], spe.sc[,2], length=0, lty=1, col='red')
#va.sc <- scores(rda_tree, choices=1:2, scaling=2, display="sp")# plots all ASVs
#text(va.sc, row.names(va.sc), cex=0.9)
#species as gray points
sp.sc <- scores(rda_tree, choices=1:2, scaling=2, display="sp")
points(sp.sc, pch=4, cex=0.5, col="gray50")
# Plot arrows of quantitative explanatory variables and their labels
spenv.sc <- scores(rda_tree, choices=1:2, scaling=2, display="bp")
arrows(0,0, spenv.sc[1:8,1], spenv.sc[1:8,2], lty=1, lwd=1.5, length=0.1, col="black")
env.names <- c("DTW", "CO2", "CH4", "DOC", "Stage_numerical")
text(spenv.sc[1:4,], env.names[1:4], cex=0.9, font=2, pos=4)
text(spenv.sc[5:8,], env.names[5:8], cex=0.9, font=2, pos=1)
#To reduce # of var:
# variance inflation factors in the RDA-- anything above 10 should be discarded
vif.cca(rda_tree) # check RDA_stepped too
# Backward selection using ordistep
str(env.sc)
RDA_stepped<-ordistep(rda(meth ~ DTW +  DOC + CO2+ m + Stage3 + CH4 +Temperature, data=env.z, direction="backward", pstep=1000, R2scop=TRUE)) #R2scope only accepts models with lower adjusted R2

summary(RDA_stepped)
# Test of RDA result
anova.cca(RDA_stepped, step=1000)
# Test of all canonical axes
anova.cca(RDA_stepped, by='axis', step=1000)
#R2, adj
(R2adj <- RsquareAdj(RDA_stepped)$adj.r.squared)# here, model explains 17% of variatio

library(ggord)
# try plotting more prettily:
ggord(RDA_stepped, metadata2$Stage) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) # looking at the raw code, this is plotting the 'wa scores', the blue dots are different species  
plot(RDA_stepped, xlab="RDA1 (15.11 %)", ylab="RDA2 (13.70 %)", 
     display=c("cn", "lc", "sp"), type="n", xlim=c(-1.3,1.3))
sites.sc <- scores(RDA_stepped, choices=1:2, scaling=2, display="lc")
text(RDA_stepped, display = "sites")

##archaeal otus for liam's permafrost data
srch<-read.csv("archaea_OTUs.csv")
methanogens<-subset(srch, Phylum=="Euryarchaeota")
colmeans(methanogens)
write.csv(methanogens,"methanogens_OTUs.csv")
