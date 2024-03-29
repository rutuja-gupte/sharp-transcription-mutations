args <- commandArgs(trailingOnly = TRUE)
# arg should be a number
fname <- args[1]
print(fname)

library(dplyr)

net.plus <- read.csv(paste0("module_data/net_plus_", fname, ".csv"))
net.minus <- read.csv(paste0("module_data/net_minus_", fname, ".csv"))
xu.tr <- read.csv("data/xu.csv")
vand.tr <- read.csv("data/vand.csv")
feat <- read.csv("data/feat.csv")

curr.chr <- 0
print("Running net.plus")
print(head(net.plus))
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
print(head(net.minus))
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

write.csv(net.plus, paste0("server_data/net_plus_", fname, ".csv"))
write.csv(net.minus, paste0("server_data/net_minus_", fname, ".csv"))

print("Job done")
