###TRY INSDICSPECIES 
library(indicspecies)
pc<-read.csv("16S data_gradient_depths.csv", row.names=1)
V4v5<-read.csv("V4V5__ASV_table.csv", row.names=1)
v6v8<-read.csv("V6V8__ASV_table.csv", row.names=1)
greenland<-read.csv("Greenland_2018__ASV_table.csv", row.names=1)
is.na(v6v8)
metaT<-read.csv("MetaT_net trap2.csv", row.names=1)
#make groupings based on site depth (numerical)
groups=c(1,rep(2,15), rep(3,6), rep(4,5))#note that samples must be arranged in order for groupings to be in order as well
#for v6v8
groups=c(1,rep(2,15), rep(3,5), rep(4,6))
#for greenland data (2 = subglacial, 3=ice marginal lake, 4= ice marginal stream)
groups=c(rep(2,2), rep(3,3), rep(4,3))
#groups for metaT
groups=c(rep(1,9), rep(2,7))
#example dataset
#data(wetland)
#groups = c(rep(1, 17), rep(2, 14), rep(3,10))
#read in metadata
meta2<-read.csv("metadata_16S_V4V5_jonesSound2019.csv", row.names=1)#doesn't work- try with groupings above, unless restructure 
#metadata to look like structured asv table
#isolate abundance data
#only if getting rid of something--abund<-V4v5[,-1]
abund<-subset(V4v5,colSums(abund)!=0)#change abund to V4v5
abund<-subset(v6v8,colSums(v6v8)!=0)
abund<-subset(greenland, colSums(greenland)!=0)
abund<-subset(metaT, colSums(metaT)!=0)
#remove any rows summing to 0
#row_sub = apply(greenland, 1, function(row) all(row !=0 ))
##Subset as usual
#a<-greenland[row_sub,]
abund<-t(abund)
head(abund)
#depth<-meta2$number_multipatt
#abund<-na.omit(abund)
str(abund)
is.na(abund)
#### now to plot
#note that duleg here makes sure that no combination of indicators (i.e, ASvs found in group 2 AND 3, and not just in either) are calculated
inv = multipatt(abund, groups, func = "r.g", duleg= T, control = how(nperm=9999)) #abund is my ASVs with species summing to 0 across the dataset removed
#groups are the sites/conditions for which I'm looking for indicators #note that index column demonstrates the group that ASV is an indicator for.
inv = multipatt(abund, groups, func = "r.g", control = how(nperm=9999)) 
#inv = multipatt(wetland, groups, func = "r.g", control = how(nperm=9999))
summary(inv)
onv2<-as.data.frame(inv$sign)
dat.multipatt.summary<-capture.output(summary(inv, indvalcomp=TRUE))
write.csv(onv2, "multipatt_depth_metaT_results.csv")
write.csv(onv2, "multipatt_depth_results_V6V8.csv")
inv3<-inv$str
##assign id to indicator species####
#read in your modified ISA table:
ISA<- read.csv("ISA_results.csv")
ISA<-read.csv("multipatt_corr_16S_indicators.csv")
#by zone
ISA<-read.csv("multipatt_zone_results.csv")
#v4v5
ISA<-read.csv("multipatt_depth_results_V4V5.csv")
#v6v8
ISA<-read.csv("multipatt_depth_results_v6v8.csv")
ISA<-as.data.table(ISA)
setkey(ISA, ASV_ID)
#read in your otu table in this way:
tax<-read.csv("Taxonomy_16S.csv")
tax<-read.csv("taxonomy_SILVA_V6V8.csv")
#RDP
tax<-read.csv("taxonomy_RDP_V4V5.csv")
tax<-as.data.table(tax)
otu2<-read.csv("V4V5__ASV_table.csv")
otu2<-as.data.table(otu2)
otu<-read.csv("16S data_gradient_depths.csv")
otu<-as.data.table(otu)
setkey(otu2, ASV_ID)
#set the key; for file named ISA, we'll be linking the ASV column together with other data files:
setkey(tax,OTU)
setkey(otu_isa, OTU)
# join the tables
otu_isa<- otu2[ISA]
#otu_isa<-otu_isa[,-c(8:29)]
#str(otu_isa$OTU)
#otu_isa$OTU<-as.factor(otu_isa$OTU)
#only if NA values are present:
#otu_isa<-na.omit(otu_isa)
#now read in taxonomy
tax<-read.csv("Taxnonomy_phyloseq.csv") #<-- change to whatever your taxonomy file is called.
#for SILVA
tax<-read.csv("V4V5_SILVA_taxonomy.csv")
tax<-read.csv("taxonomy_RDP_V4V5.csv")
tax<-read.csv("taxonomy_RDP_V6V8.csv")
otu<-read.csv("V4V5__ASV_table.csv")
otu<-read.csv("V6V8__ASV_table.csv")
otu<-as.data.table(otu)
tax<-as.data.table(tax)
setkey(tax, ASV_ID)
otu_isa2 <- tax[otu_isa]
#remove col 10:12
otu_isa2<-otu_isa2[,-c(7:8)]
otu_isa2<-na.omit(otu_isa2)
write.csv(otu_isa2,"multipatt_ISA_corr.csv")
write.csv(otu_isa2, "multipatt_ISA_16S_tax_otu_RDP.csv")
write.csv(otu_isa2, "mulitpatt_ISA_tax_otu_RDP_V6V8.csv")
#write this csv file to your directory, and edit according to read me.
write.csv(otu_isa2, "otu_taxa_isa.csv")
otu_isa2
#### plotting this data #######
#v4v5
isa_silva_v4v5<-read.csv("multipatt_ISA_16S_tax_otu_SILVA_V4V5.csv")
#RDP
isa_RDP_v4v5<-read.csv("multipatt_ISA_16S_tax_otu_RDP_V4V5.csv")
#make long
long_v4v5<-gather(isa_silva_v4v5, site, count, Ter_122:VIO_10_95m)
#RDP
long_rdp<-gather(isa_RDP_v4v5, site, count, Ter_122:VIO_10_95m)
#aggregate:
ag_v4v5<-aggregate(count~Order +site+index,data=long_rdp,FUN=sum)
ag_v4v5<-as.data.table(ag_v4v5)
#make wide:
wide_v4v5<-reshape(ag_v4v5, idvar="Order", timevar="site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide_v4v5)){
  colnames(wide_v4v5)[col] <-  sub("Count.", "", colnames(wide_v4v5)[col])
}
wide_v4v5<-as.data.table(wide_v4v5)
setkey(wide_v4v5, Order)
#v6v8
isa_silva_v6v8<-read.csv("mulitpatt_ISA_tax_otu_RDP_V6V8.csv")
#RDP
#make long
long_v6v8<-gather(isa_silva_v6v8, site, count, o_122:X48.o)
#aggregate:
ag_v6v8<-aggregate(count~Order +site+index,data=long_v6v8,FUN=sum)
#write,
write.csv(ag_v6v8, "V6V8_RDP_aggregated_Indicator_orders.csv")
#make wide:
wide_v6v8<-reshape(ag_v6v8, idvar="Order", timevar="site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide_v6v8)){
  colnames(wide_v6v8)[col] <-  sub("Count.", "", colnames(wide_v6v8)[col])
}
ag_v6v8<-as.data.table(ag_v6v8)
setkey(ag_v6v8, Order)
#merge v4v5 & v6v8
total_isa<-wide_v4v5[wide_v6v8]
#export
write.csv(total_isa, "combined_indicator_Order_count_RDP_v8v6_v4v5.csv")
#read in
total_isa<-read.csv("combined_indicator_Order_count_RDP_v8v6_v4v5.csv")
#For greenland data#----------------
greenland<-read.csv("multipatt_type_results_significant_ind.csv")

#make long
long<-gather(greenland,Site, Rel, B:L )
###
long<-gather(total_isa,site, count, Ter_122:X81.o )
write.csv(long, "long_format_combined_indicators_RDP_v6v8v4v5.csv")
#after organizing in excel:
new<-read.csv("long_format_combined_indicators_SILVA_v6v8v4v5.csv")
#rdp
ewn<-read.csv("long_format_combined_indicators_RDP_v6v8v4v5.csv")
#remove zero counts:
new2<-subset(new, count>0)
new2<-subset(ewn, count>0)
#greenland data#
new<-subset(long, Rel>0)
#calc rel abund:
new2$rel_abund=(new2$count/sum(new2$count))
#find range of rel_abund
range(new2$rel_abund)
range(new2$count)
#filter out, for RDP certain orders in exce;:
new2$sqrt_count<-sqrt(new2$count)
#sweet- now plot (note that this is for raw counts; change for relative abundance)
theme_set(theme_bw())
col=c( "subglacial stream"="darkorchid3", "50m<x>100m"="goldenrod")## <- change accordingly, if desired-- see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
#order silva:
new2$Order<-factor(new2$Order, levels=c("Acetobacterales", "Babeliales","Bacilalles", "Betaproteobacteriales", "Candidatus Kaiserbacteria", "Candidatus Nomurabacteria",
                                        "Candidatus Pacebacteria", "Chitinophagales", "Chthoniobacterales", "Diplorickettsiales", "Fibrobacterales", "Frankiales", "Gaiellales",
                                        "Leptolyngbyales", "Micrococcales", "Microtrichales", "Ochrophyta", "Paracaedibacterales", "Pedosphaerales", "Saccharimonadales", "Solibacterales",
                                        "Tepidisphaerales", "uncultured Halanaerobiaceae bacterium", "Bacillales", "Bacteroidales", "Bdellovibrionales", "Caulobacterales", "Chlamydiales",
                                        "Chloroplast", "Clostridiales", "Cytophagales", "Desulfuromonadales", "Flavobacteriales", "Gemmatimonadales", "Ignavibacteriales",
                                        "Myxococcales", "Oceanospirillales", "Oligoflexales", "Pseudomonadales", "Rhizobiales", "Rhodobacterales", "Rickettsiales", "Sphingobacteriales", "Sphingomonadales",
                                        "Xanthomonadales"))





#order RDP:
new2$Order<-factor(new2$Order, levels=c("Acidimicrobidae", "Actinobacteridae", "Aerolineales", "Aridibacter", "Bryobacter", "Burkholderiales", "Chthonomodales", "Coriobacteridae",
                                        "Deinococcales", "Elusimicrobiales", "Gallionellales", "Granulicella", "Holophagales", "Ktedonobacterales", "Legionellales", "Methylophilales", "Neisseriales", "Nitrosomodales",
                                        "Nitrospirales", "Opitutales", "Paludibaculum", "Phycisphaerales", "Planctomycetales", "Puniceicoccales", "Rhodospirillales","Rubrobacteridae",
                                        "Selenomodales", "Sulfuricellales", "Victivallales", "Bacillales", "Bacteroidales", "Bdellovibrioles", "Caulobacterales", "Chlamydiales", "Chloroplast", "Clostridiales",
                                        "Cytophagales", "Desulfuromodales", "Flavobacteriales", "Gemmatimodales", "Igvibacteriales", "Myxococcales", "Oceanospirillales", "Oligoflexales", "Pseudomodales",
                                        "Rhizobiales","Rhodobacterales", "Rickettsiales","Sphingobacteriales", "Sphingomodales", "Xanthomodales"))

#remove NA
new2<-na.omit(new2)
#for greenland data

new$Index<-factor(new$Index, levels=c("subglacial", "ice_marginal_stream", "ice_marginal_lake"))
col=c("subglacial"="grey", "ice_marginal_lake"="lightblue", "ice_marginal_stream"="tan")

g <- ggplot(new, aes(x=Site, y=Family, size=Rel, colour=Index)) 
g + geom_point() +
  scale_colour_manual(values=col)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10))+
  theme(axis.text.y = element_text( size=9))+
  theme(axis.title.x=element_text(size=12))+
  theme(axis.title.y=element_text(size=12))+ggtitle('Indicator organisms')+
  scale_fill_manual(name="relative abundance")
  
#for outlining points in black https://t-redactyl.io/blog/2016/02/creating-plots-in-r-using-ggplot2-part-6-weighted-scatterplots.html



















#for other projects, how I plotted:
#for straight up plotting raw counts, use:
otu_type<-read.csv("for plot_isa.csv") # make sure its the modified version- see readme
#For plotting relative abundance:
otu_rel<-read.csv("rel_abund_isa.csv", row.names=1) # make sure its the modified version- see readme
totu_type<-t(otu_rel)
x<-long/rowSums(totu_type)
x<-na.omit(x)
x<-t(x)
x<-as.data.frame(x)
write.csv(x, "rel-abund_isa2.csv")
#read this file back in once you've edited it according to readme.
rel_abund_isa<-read.csv("rel-abund_isa2.csv")
#Make your raw counts table a long format table:
x_long<-gather(otu_type, Site, Count,MB135_J:YB45_S)# make sure you align this to your first sample name: to your last , i.e., MB135_J:YB45_S )
# do the same for your  relative abundance table
rel_long<-gather(rel_abund_isa, Site, Rel.Abund,MB135:YB10)# make sure you align this to your first sample name: to your last , i.e., MB135_J:YB45_S )
#write both as csv files into your directory, ex:
write.csv(x_long, "for_plot_rel_abund.csv")
#remove 0s
for_plot <- filter(x_long, Count!= 0)
#for relative abundance:
for_plot_rel <- filter(rel_long, Rel.Abund!= 0)
#order according to group level (do the same for relative abundance table, changing the appropriate file names)
for_plot$name<-as.factor(for_plot$name)
#order
for_plot$name<-factor(for_plot$name, levels=c( "river", "natural biofilm", "artificial biofilm"))
#Now order your x axis (site names) <-- edit accordingly.
for_plot$Site <- factor(for_plot$Site, levels = c( "P510_J","P2030_J","YB05_J", "YB510_J", "YB2030_J", "YB3040_J", "YB4050_J", "YB125135_J",
                                                   "YB150160_J", "MB05_J", "MB2030_J", "MB7080_J", "MB8090_J", "MB125135_J", "MB150160_J"))


#for relative abundance x axis ordering:
for_plot_rel$Site <- factor(for_plot_rel$Site, levels = c("YB5_J", "YB10_J", "YB25_J", "YB35_J", "YB45_J", "YB135_J",
                                                          "YB155_J", "MB5_J", "MB15_J", "MB25_J","MB45_J", "MB75_J", "MB85_J", "MB135_J", "MB155_J","YB5_S", "YB10_S", "YB25_S", "YB35_S", "YB45_S", "YB135_S",
                                                          "YB155_S", "MB5_S","MB15_S", "MB25_S","MB45_S", "MB55_S", "MB75_S", "MB85_S", "MB135_S", "MB155_S"))
# read

#sweet- now plot (note that this is for raw counts; change for relative abundance)
theme_set(theme_bw())

col=c( "chartreuse4","darkgoldenrod1", "blue")## <- change accordingly, if desired-- see http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
g <- ggplot(for_plot, aes(x=Site, y=Family, colour=group)) 
g + geom_point(aes(size=Count)) +
  scale_colour_manual(values=col)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9))+
  theme(axis.text.y = element_text( size=9))+
  theme(axis.title.x=element_text(size=12))+
  theme(axis.title.y=element_text(size=12))
#export both bubble plots as PDF.