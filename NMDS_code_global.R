##NMDS
#Steps
#1. first must aggregate all common orders or families, across all sites
#2. do for other dataset
#3. then make datasets same size
#4. then merge based on shared orders and run NMDS, colouring based on primer
##### accomplishing steps #####
#1: based on SILVA: 
v4<-read.csv("V4V5__ASV_table.csv")
v4<-read.csv("V4V5_ASV_table_for_phytoref.csv")
#---------------------------------#other example data---------------------
#perma<-read.csv("Original_non_transformed_OTU_table_global.csv")
#NSR_2019<-read.csv("NSR_ASV_table_raw_2019.csv")
#NSR_2020<-read.csv("NSR_2020_ASV_table.csv")
#CMN<-read.csv("16S_CMN2019_table_Maria_working_copy.csv")
#js<-read.csv("16S_JonesSound2019_ASV.csv")
#ext_test<-read.csv("Extraction_test_ASV.csv")
#----------------------------------------------------------------------
#remove first column
v4<-v4[,-1]
#---------------------------------#other example data---------------------
#perma<-perma[,-1]
#nsr_2019<-NSR_2019[,-1]
#nsr_2020<-NSR_2020[,-1]
#cmn<-CMN[,-1]
#js<-js[,-1]
#ext_test<-ext_test[,-1]
#----------------------------------------------------------------------
#make long
long_v4<-gather(v4, Site, Count, Ter_122:VIO_10_95m)
#---------------------------------#other example data---------------------
#perma_long<-gather(perma, Site, Count,OB_125135_Jun:YS_4050)
#nsr2019_long<-gather(nsr_2019, Site, Count, X1GM:White)
#nsr2020_long<-gather(nsr_2020, Site, Count, X1:Q3)
#cmn_long<-gather(cmn, Site, Count, X4:SRLO.16S)
#js_long<-gather(js, Site, Count, VIO_27_3m:VIO_31_20m)
#ext_test_long<-gather(ext_test, Type, Count, M:L)
#----------------------------------------------------------------------
#number of unique orders:
length(unique(long_v4$class))
#merge all duplicate orders by summing up counts belonging to orders: 
new<-aggregate(Count~class +Site,data=long_v4,FUN=sum)
#turn new into wide format:
wide<-reshape(new, idvar="class", timevar="Site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide)){
  colnames(wide)[col] <-  sub("Count.", "", colnames(wide)[col])
}
wide<-wide[,-1] # you may not need to do this...
#write.csv
write.csv(wide, "V4V5_Phytoref_classes_ASV.csv")
#---------------------------------#other example data---------------------
#write.csv(wide, "permafrost_global_ASV.csv")
#write.csv(wide, "NSR_2019_global_ASV.csv")
#write.csv(wide, "NSR_2020_global_ASV.csv")
#write.csv(wide, "cmn_global_ASV.csv")
#write.csv(wide, "js_2019_global_ASV.csv")
#write.csv(wide, "ext_test_global_Orders_greengenes.csv")
#----------------------------------------------------------------------
#Do for remaining datasets, i.e., here for the v6v8 dataset:
v6<-read.csv("V6V8__ASV_table.csv")
v6<-read.csv("V6V8_ASV_table_phytoref.csv")
#for knorr project:
taxa<-read.csv("16S data.csv")
#remove otu ID and other taxa columns that aren't  phylum
taxa<-taxa[-c(1,2,4,5,6,7,8),]
#make long
long<-gather(taxa, Site, Count, KN_S15_1000m:KN_S7_NADW)
#merge all duplicate orders by summing up counts belonging to orders: 
new<-aggregate(Count~Phylum +Site,data=long,FUN=sum)
#turn new into wide format:
wide<-reshape(new, idvar="Phylum", timevar="Site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide)){
  colnames(wide)[col] <-  sub("Count.", "", colnames(wide)[col])
}
#write.csv
write.csv(wide, "Phyla_knorr_16S.csv")
##### for Orion's dataset:#################
#remove first column
v6<-v6[,-1]
#make long
long_v6<-gather(v6, Site, Count, Ter_122:VIO_10_95m)

#number of unique orders:
length(unique(long_v6$Class))
#merge all duplicate orders by summing up counts belonging to orders: 
new<-aggregate(Count~Class +Site,data=long_v6,FUN=sum)
#turn new into wide format:
wide<-reshape(new, idvar="Class", timevar="Site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide)){
  colnames(wide)[col] <-  sub("Count.", "", colnames(wide)[col])
}
#write.csv
write.csv(wide, "V6V8_Order_Count_RDP.csv")
write.csv(wide,"V6V8_Phytoref_classes_ASV.csv")
#try to merge both new datasets by shared order:
v4<-read.csv("V4V5_Order_Count_RDP.csv")
v4<-as.data.table(v4)
setkey(v4, Order)
v6<-read.csv("V6V8_Order_Count_RDP.csv")
v6<-as.data.table(v6)
setkey(v6, Order)
#merge
v4_v6<-v4[v6]
#merge multiple datasets:
multmerge = function(path){
  filenames=list.files(path=path, full.names=TRUE)
  rbindlist(lapply(filenames, fread))
}#doesn't work- must have even columns across all datasets
path <- "Dropbox (Bhatia Lab)/Bhatia Lab Team Folder/Projects/Devon/Devon_2019/Data/16S/R_analysis/global_NMDS_files"
DF <- multmerge(path)
# you may need to do this manually as per steps below:
perma<-read.csv("permafrost_global_ASV.csv")
perma<-perma[,-1]
perma<-as.data.table(perma)
setkey(perma, Order)
cmn<-read.csv("cmn_global_ASV.csv")
cmn<-cmn[,-1]
cmn<-as.data.table(cmn)
setkey(cmn, Order)
#merge
perma_cmn<-perma[cmn]
setkey(perma_cmn, Order)
nsr<-read.csv("NSR_2019_global_ASV.csv")
nsr<-nsr[,-1]
nsr<-as.data.table(nsr)
setkey(nsr, Order)
#merge
perma_cmn_nsr2019<-nsr[perma_cmn]
setkey(perma_cmn_nsr2019, Order)
nsr2<-read.csv("NSR_2020_global_ASV.csv")
nsr2<-nsr2[,-1]
nsr2<-as.data.table(nsr2)
setkey(nsr2, Order)
#merge
perma_cmn_nsr20192020<-nsr2[perma_cmn_nsr2019]
setkey(perma_cmn_nsr20192020, Order)
js<-read.csv("js_2019_global_ASV.csv")
js<-js[,-1]
js<-as.data.table(js)
setkey(js, Order)
perma_cmn_nsr20192020_js<-js[perma_cmn_nsr2019]
#merge
#remove orders not shared:
all<-na.omit(perma_cmn_nsr20192020_js)
#only for SILVA dataset: all<-all[,-1]
write.csv(all, "combined_global_SILVA_ASV.csv")
all<-read.csv("Book2.csv", row.names=1)
all<-t(all)
meta<-read.csv("metadata_global.csv", row.names=1)
#################for an nmds-------------------------------------------------------------------
#get grouping info:
grouping_info<-data.frame(row.names=rownames(all),t(as.data.frame(strsplit(rownames(all),"X"))))
#calculate NMDS matrix
sol<-metaMDS(all,distance = "bray", k = 2, trymax = 50)
#make dataframe:

#plot

col=c("blue","red")
p<-ggplot(data=NMDS,aes(x,y,colour=depth))+
  geom_point(size=8)+
  scale_colour_manual(values=c(col))+
  geom_point(shape=1, size=8, colour="black")
p + theme_bw()
#significance
ano=anosim(all, NMDS$depth, distance="bray", permutations=999)

NMDS=data.frame(x=sol$point[,1],y=sol$point[,2], Site_v6v8=as.factor(grouping_info[,1]),Site_v4v5=as.factor(grouping_info[,2]))
NMDS$primer<-c("V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5",
               "V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V4V5","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8",
               "V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8","V6V8")
col=c("blue","pink")

p<-ggplot(data=NMDS,aes(x,y,colour=primer))+
  geom_point(aes(colour=primer), size=8)+
  geom_point(shape=1, size=8, colour='black')+
  scale_colour_manual(values=col)+
  labs(title= "16S NMDS, based on RDP-assigned Orders", colour= "primer")+
  xlab("NMDS1")+
  ylab("NMDS2")
p + theme_bw()
range(NMDS$x)