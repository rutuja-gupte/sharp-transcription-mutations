library(dplyr)
library(profvis)

net.plus <- read.csv("data/test_plus.csv")[1:40000,]
net.minus <- read.csv("data/test_minus.csv")[1:40000,]
xu.tr <- read.csv("data/xu.csv")
vand.tr <- read.csv("data/vand.csv")
feat <- read.csv("data/feat.csv")

merged = bind_rows(xu.tr,vand.tr,feat)
merged.plus = merged[merged$strand == "+",]
merged.minus = merged[merged$strand == "-",]

net.plus.pos.chr = paste(net.plus$pos,net.plus$chr_num,sep="-")
net.plus.minus.chr = paste(net.minus$pos,net.minus$chr_num,sep="-")

# profvis({

search.sign = function(pos.chr,left,right,strand,chr_num){
  split = as.numeric(strsplit(pos.chr,"-")[[1]])
  pos = split[1]
  chr = split[2]
  
  chr_num == chr & left <= pos & right >= pos
}


curr.chr <- 0
print("Running net.plus")
for (i in 1:nrow(net.plus)){
  chr <- net.plus[i,"chr_num"]
  if (curr.chr != chr) {
    print(Sys.time())
    print("Now on:")
    print(chr)
    curr.chr = chr
  }
  pos <- net.plus[i,"pos"]
  a <- xu.tr[xu.tr$chr_num == chr & xu.tr$strand == "+" & 
               xu.tr$left <= pos & xu.tr$right >= pos,]
  b <- vand.tr[vand.tr$chr_num == chr & vand.tr$strand == "+" & 
                 vand.tr$left <= pos & vand.tr$right >= pos,]
  c <- feat[feat$chr_num == chr & feat$strand == "+" & 
              feat$left <= pos & feat$right >= pos,]
  
  df <- bind_rows(a,b,c)
  if (nrow(df) != 0){
    df$diff <- pos - df$left
    mindiff <- min(df$diff)
    vals <- df[df$diff == mindiff,][1,]
    net.plus[i,]$net_tot <- vals$net_tot
    net.plus[i,]$net_dens <- vals$net_dens
    net.plus[i,]$class <- vals$class
  }
}

curr.chr <- 0
print("Running net.minus")
for (i in 1:nrow(net.minus)){
  chr <- net.minus[i,"chr_num"]
  if (curr.chr != chr) {
    print(Sys.time())
    print("Now on:")
    print(chr)
    curr.chr = chr
  }
  pos <- net.minus[i,"pos"]
  a <- xu.tr[xu.tr$chr_num == chr & xu.tr$strand == "-" & 
               xu.tr$left <= pos & xu.tr$right >= pos,]
  b <- vand.tr[vand.tr$chr_num == chr & vand.tr$strand == "-" & 
                 vand.tr$left <= pos & vand.tr$right >= pos,]
  c <- feat[feat$chr_num == chr & feat$strand == "-" & 
              feat$left <= pos & feat$right >= pos,]
  
  df <- bind_rows(a,b,c)
  if (nrow(df) != 0){
    df$diff <- pos - df$left
    mindiff <- min(df$diff)
    vals <- df[df$diff == mindiff,][1,]
    net.minus[i,]$net_tot <- vals$net_tot
    net.minus[i,]$net_dens <- vals$net_dens
    net.minus[i,]$class <- vals$class
  }
}

write.csv(net.plus, "test_plus_sites2.csv")
write.csv(net.minus, "test_minus_sites2.csv")


# })


