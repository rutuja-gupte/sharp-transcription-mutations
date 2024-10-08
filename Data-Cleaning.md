Data Processing
================
Rutuja
2024-09-04

This document contains some code that was originally written by Dr
Nathaniel Sharp.

``` r
# Import
library(tidyverse)
```

``` r
# Key for replacing chromosome names later

chr_names<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

chr_names2<-c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16")
```

## Data

### NET seq dataset

**Source:** Churchman, L. Stirling, and Jonathan S. Weissman. “Nascent
Transcript Sequencing Visualizes Transcription at Nucleotide
Resolution.” Nature 469, no. 7330 (January 2011): 368–73.
<https://doi.org/10.1038/nature09652>.

The data comes from deep sequencing 40 bases from the 3’ ends of the
nascent transcripts isolated from the cells. The dataset has single
nucleotide resolution.

**Strain:** BY4741, mid exponential phase  
genotype/variation: RPB3-3xFLAG::NAT his3D1::TEF2-GFP-Adh1 KAN  
library strategy: NET-Seq  
molecule: nascent RNA  
ip: IP against 3xFLAG affinity label”

Here we are reading in the plus and minus datasets, adding chromosome
numbers. Then converting the number of reads back to integers. Finally,
combining the datasets into one table with strands labelled.

``` r
net.plus<-read.delim("transcription_data/GSM617027_WT_NC_plus.txt",h=F,sep=" ")
names(net.plus)<-c("pos","value")
prev<-c(0,net.plus$pos[1:(nrow(net.plus)-1)])
diff<-sign(net.plus$pos-prev)
net.plus$chr<-NA
cut.rows<-c(1,which(diff<0),nrow(net.plus))
for(i in c(1:16)){net.plus$chr[cut.rows[i]:cut.rows[i+1]]<-i}


net.minus<-read.delim("transcription_data/GSM617027_WT_NC_minus.txt",h=F,sep=" ")
names(net.minus)<-c("pos","value")
prev<-c(0,net.minus$pos[1:(nrow(net.minus)-1)])
diff<-sign(net.minus$pos-prev)
net.minus$chr<-NA
cut.rows<-c(1,which(diff<0),nrow(net.minus))
for(i in c(1:16)){net.minus$chr[cut.rows[i]:cut.rows[i+1]]<-i}

rm(diff, prev, cut.rows)
```

``` r
val.plus <- min(net.plus$value)
val.minus <- min(net.minus$value)

# re-scale y so that it is an integer (assumes lowest non-zero value of y corresponds to 1 read per 10^7 sequences)
net.plus$reads<-net.plus$value/val.plus
net.minus$reads<-net.minus$value/val.minus
```

``` r
net.plus$strand <- "+"
net.minus$strand <- "-"

net <- bind_rows(net.plus, net.minus)
rm(net.plus, net.minus)
```

### Power to detect mutations dataset

The dataset is based on coverage thresholds from sequencing of MA lines.
Includes the reference genome. I am excluding regions with 0 power.
Finally combining information with the netseq data. We probably do not
need the power later but it is good to have that information. I mainly
need the reference genome later.

``` r
# power to detect mutations
power <- read.delim("power_by_site_with_bases.txt", h=F)
names(power) <- c("chr", "pos", "hap.wt", "hap.del", "dip.wt", "dip.del", "ref")
power$chr<- sapply(power$chr, function(i){which(chr_names2 == i)})
power <- power[!(power$hap.wt == 0 & power$hap.del == 0 &
                   power$dip.wt == 0 & power$dip.del == 0), ]
power1 <- power
power1$strand <- "+"
power2 <- power
power2$strand <- "-"
power <- bind_rows(power1, power2)
rm(power1, power2)
```

``` r
d <- left_join(power, net, by=c("chr", "pos", "strand")) %>% replace_na(list("value"=0, "reads"=0))
rm(net, power)
```

Observation: transcription data from about 20% of the genome

### Xu dataset

``` r
# xu.headers <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "sourceMeta")
# 
# xu.ofr <- read.delim("xu_data/Xu_2009_ORF-Ts_V64.txt", header=F)
# names(xu.ofr) <- xu.headers
# 
# xu.cuts <- read.delim("xu_data/Xu_2009_CUTs_V64.txt", header=F)
# names(xu.cuts) <- xu.headers
# 
# xu.suts <- read.delim("xu_data/Xu_2009_SUTs_V64.txt", header=F)
# names(xu.suts) <- xu.headers
# 
# xu.others <- read.delim("xu_data/Xu_2009_other_transcripts_V64.txt", header=F)
# names(xu.others) <- xu.headers
# 
# xu <- bind_rows(xu.ofr, xu.cuts, xu.suts, xu.others)
# rm(xu.ofr, xu.cuts, xu.suts, xu.others)
# 
# xu$left<-xu$start
# xu$right<-xu$end
# xu$chr_num<-sapply(xu$chr,function(x){which(chr_names==x)})
# xu$class <- data.frame(str_split_fixed(xu$source, "_", 3))$X3
# xu <- xu %>% select(chr_num, left, right, strand, class)
# xu.tr <- xu
# rm(xu)
# 
# vand.tr<-read.delim("transcription_data/vanD_annotations.txt",h=F)
# vand.tr<-vand.tr[,c(1,4,5,7)]
# names(vand.tr)<-c("chr_name","left","right","strand")
# vand.tr$chr_num<-sapply(vand.tr$chr_name,function(i){which(chr_names==i)})
# vand.tr <- vand.tr %>% select(-chr_name)
# 
# feat<-read.csv("transcription_data/features.csv",h=T)
# feat<-feat[,c(2:5,10:13)]
# names(feat)<-c("SGDID","feat_type","feat_qual","name", "chr_num", "start", "stop", "strand_name")
# feat<-feat[feat$feat_qual != "Dubious",]
# 
# feat$chr_num<-as.numeric(feat$chr_num)
# feat <- feat %>% filter(!is.na(chr_num))
# 
# feat<-feat[feat$chr_num %in% c(1:16),]
# feat$strand<-sapply(feat$strand_name,function(i){if(i=="W"){"+"}
#   else{
#     if(i=="C"){"-"}
#     else{NA}}})
# feat <- feat %>% filter(!is.na(strand))
# feat$left<-sapply(1:nrow(feat),function(r){if(feat$strand[r]=="+"){feat$start[r]}else{feat$stop[r]}})
# feat$right<-sapply(1:nrow(feat),function(r){if(feat$strand[r]=="+"){feat$stop[r]}else{feat$start[r]}})
# 
# # needs to be modified if I missed something that gets transcribed or if I added something that does not get transcribed.
# feat_types <- c("ORF", "intron", "ncRNA_gene", "noncoding_exon", "tRNA_gene", "snoRNA_gene", "pseudogene", "five_prime_UTR_intron", "snRA_gene", "rRNA_gene", "telomerase_RNA_gene", "CDS", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "transposable_element_gene")
# 
# feat <- feat %>% filter(feat_type %in% feat_types)
# feat <- feat %>% select(feat_type, chr_num, left, right, strand)
# 
# d <- power %>% select(chr_num, pos, ref, hap.wt, hap.del, dip.wt, dip.del)
# d$pwr.tot<-rowSums(d[,c(4:7)])
# 
# xu.tr <- xu.tr %>% relocate(chr_num, left, right, strand)
# vand.tr$class <- "XUTs"
# vand.tr <- vand.tr %>% relocate(chr_num, left, right, strand)
# feat <- feat %>% rename("class" = "feat_type")
# feat <- feat %>% relocate(chr_num, left, right, strand)
# 
# feats <- bind_rows(xu.tr, vand.tr, feat)
# rm(xu.tr, vand.tr, feat)
# rm(power)
# 
# net.minus$strand <- "-"
# net.plus$strand <- "+"
# net <- bind_rows(net.plus, net.minus)
# rm(net.plus, net.minus)
# 
# # now write the clean things to clean data files
# write.csv(d, "clean_data/d.csv")
# write.csv(feats, "clean_data/features.csv")
# write.csv(net, "clean_data/net.csv")
```
