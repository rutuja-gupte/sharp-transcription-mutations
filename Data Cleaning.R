library(tidyverse)

chr_names<-c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII", "chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

chr_names2<-c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12","chr13","chr14","chr15","chr16")

net.plus<-read.delim("transcription_data/GSM617027_WT_NC_plus.txt",h=F,sep=" ")
names(net.plus)<-c("pos","value")
net.minus<-read.delim("transcription_data/GSM617027_WT_NC_minus.txt",h=F,sep=" ")
names(net.minus)<-c("pos","value")

prev<-c(0,net.plus$pos[1:(nrow(net.plus)-1)])
diff<-sign(net.plus$pos-prev)
net.plus$chr_num<-NA
cut.rows<-c(1,which(diff<0),nrow(net.plus))
for(i in c(1:16)){net.plus$chr_num[cut.rows[i]:cut.rows[i+1]]<-i}

prev<-c(0,net.minus$pos[1:(nrow(net.minus)-1)])
diff<-sign(net.minus$pos-prev)
net.minus$chr_num<-NA
cut.rows<-c(1,which(diff<0),nrow(net.minus))
for(i in c(1:16)){net.minus$chr_num[cut.rows[i]:cut.rows[i+1]]<-i}

rm(diff, prev, cut.rows)

val.plus <- min(net.plus$value)
val.minus <- min(net.minus$value)

# re-scale y so that it is an integer (assumes lowest non-zero value of y corresponds to 1 read per 10^7 sequences)
net.plus$reads<-net.plus$value/val.plus
net.minus$reads<-net.minus$value/val.minus

# power to detect mutations
power<-read.delim("power_by_site_with_bases.txt",h=F)
names(power)<-c("chr","pos","hap.wt","hap.del","dip.wt","dip.del","ref")
power$chr_num<-sapply(power$chr,function(i){which(chr_names2==i)})

xu.headers <- c("chr", "source", "type", "start", "end", "score", "strand", "frame", "sourceMeta")

xu.ofr <- read.delim("xu_data/Xu_2009_ORF-Ts_V64.txt", header=F)
names(xu.ofr) <- xu.headers

xu.cuts <- read.delim("xu_data/Xu_2009_CUTs_V64.txt", header=F)
names(xu.cuts) <- xu.headers

xu.suts <- read.delim("xu_data/Xu_2009_SUTs_V64.txt", header=F)
names(xu.suts) <- xu.headers

xu.others <- read.delim("xu_data/Xu_2009_other_transcripts_V64.txt", header=F)
names(xu.others) <- xu.headers

xu <- bind_rows(xu.ofr, xu.cuts, xu.suts, xu.others)
rm(xu.ofr, xu.cuts, xu.suts, xu.others)

xu$left<-xu$start
xu$right<-xu$end
xu$chr_num<-sapply(xu$chr,function(x){which(chr_names==x)})
xu$class <- data.frame(str_split_fixed(xu$source, "_", 3))$X3
xu <- xu %>% select(chr_num, left, right, strand, class)
xu.tr <- xu
rm(xu)

vand.tr<-read.delim("transcription_data/vanD_annotations.txt",h=F)
vand.tr<-vand.tr[,c(1,4,5,7)]
names(vand.tr)<-c("chr_name","left","right","strand")
vand.tr$chr_num<-sapply(vand.tr$chr_name,function(i){which(chr_names==i)})
vand.tr <- vand.tr %>% select(-chr_name)

feat<-read.csv("transcription_data/features.csv",h=T)
feat<-feat[,c(2:5,10:13)]
names(feat)<-c("SGDID","feat_type","feat_qual","name", "chr_num", "start", "stop", "strand_name")
feat<-feat[feat$feat_qual != "Dubious",]

feat$chr_num<-as.numeric(feat$chr_num)
feat <- feat %>% filter(!is.na(chr_num))

feat<-feat[feat$chr_num %in% c(1:16),]
feat$strand<-sapply(feat$strand_name,function(i){if(i=="W"){"+"}
  else{
    if(i=="C"){"-"}
    else{NA}}})
feat <- feat %>% filter(!is.na(strand))
feat$left<-sapply(1:nrow(feat),function(r){if(feat$strand[r]=="+"){feat$start[r]}else{feat$stop[r]}})
feat$right<-sapply(1:nrow(feat),function(r){if(feat$strand[r]=="+"){feat$stop[r]}else{feat$start[r]}})

# needs to be modified if I missed something that gets transcribed or if I added something that does not get transcribed.
feat_types <- c("ORF", "intron", "ncRNA_gene", "noncoding_exon", "tRNA_gene", "snoRNA_gene", "pseudogene", "five_prime_UTR_intron", "snRA_gene", "rRNA_gene", "telomerase_RNA_gene", "CDS", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "transposable_element_gene")

feat <- feat %>% filter(feat_type %in% feat_types)
feat <- feat %>% select(feat_type, chr_num, left, right, strand)

d <- power %>% select(chr_num, pos, ref, hap.wt, hap.del, dip.wt, dip.del)
d$pwr.tot<-rowSums(d[,c(4:7)])

xu.tr <- xu.tr %>% relocate(chr_num, left, right, strand)
vand.tr$class <- "XUTs"
vand.tr <- vand.tr %>% relocate(chr_num, left, right, strand)
feat <- feat %>% rename("class" = "feat_type")
feat <- feat %>% relocate(chr_num, left, right, strand)

feats <- bind_rows(xu.tr, vand.tr, feat)
rm(xu.tr, vand.tr, feat)
rm(power)

net.minus$strand <- "-"
net.plus$strand <- "+"
net <- bind_rows(net.plus, net.minus)
rm(net.plus, net.minus)

# now write the clean things to clean data files
write.csv(d, "clean_data/d.csv")
write.csv(feats, "clean_data/features.csv")
write.csv(net, "clean_data/net.csv")
