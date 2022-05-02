
library(stringr)

NDD_genes <- read.csv("Data/Repeats_of_interest_NDD_ID.txt", sep="\t", header=T)
Other_genes <- read.csv("Data/Repeats_of_interest_Other.txt", sep="\t", header=T)
No_genes <- read.csv("Data/Repeats_of_interest_No_association.txt", sep="\t", header=T)

ID_panel <- read.csv("C:/Users/DLE/Dropbox/University/PhD/Writing/Manuscript 1/CGG Catalogue Paper/Misc/ID_Genes_list.txt", header=F)$V1

NDD_genes$In_panel <- apply(NDD_genes, 1, 
                            function(x){
    spl <- str_split(x[1], "->")
    if (spl[[1]][1] %in% ID_panel | spl[[1]][1] %in% ID_panel){T} else{F}
})
Other_genes$In_panel <- apply(Other_genes, 1, 
                            function(x){
                              spl <- str_split(x[1], "->")
                              if (spl[[1]][1] %in% ID_panel | spl[[1]][1] %in% ID_panel){T} else{F}
                            })
No_genes$In_panel <- apply(No_genes, 1, 
                            function(x){
                              spl <- str_split(x[1], "->")
                              if (spl[[1]][1] %in% ID_panel | spl[[1]][1] %in% ID_panel){T} else{F}
                            })

test2 <- merge(NDD_genes[c("Gene", "Region")], tplot, by=c("Gene", "Region"))

