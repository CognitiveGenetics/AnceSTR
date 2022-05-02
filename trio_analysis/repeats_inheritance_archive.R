#     INTERGENERATIONAL REPEATS ANALYSIS
#     IMPORT FUNCTIONS
source("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/anceSTR/anceSTR/functions/functions.R")
#     LIBRARIES
Package_up()
#     IMPORT BASE DATA
c(ngc, loci) %<-% Import_Data(as.character(Sys.info()[['sysname']]), 38)
ngc <- read.csv(tcltk::tk_choose.files(caption = "Select Reads File"), header = T, sep = "\t")
loci <- read.csv(tcltk::tk_choose.files(caption = "Select Corresonding Loci File"), header = T, sep = "\t")
#ngc <- test("D:/Work_Storage/MSSNG/Trios/CGG_Repeats_Trios_20220215_filtered.csv","C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Sample Data/MSSNG_family_manifest.tab", 3)
ngc <- assign_structure("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Analysis/Quads/CGG_Repeats_Quads_Full_20220215_filtered.csv","C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Sample Data/MSSNG_family_manifest.tab", 4)
ngc <- read.csv("D:/Work_Storage/MSSNG/Trios_MSSNG/CGG_Repeats_Full_Trios_filtered.tab", header = T, sep = "\t")

#     ORGANISE RAW DATA INTO TRIOS
c(trios_df, families_df) %<-% Get_Trios(ngc)
#     ORGANISE AND CLASSIFY TRIO READS
c(Events_df, Poly_df, Stable_df, mdata) %<-% Collect_Expansion_Events(families_df, "trios")
#     ANALYSE REPEAT LENGTH CHANGES
#     Determine which are the Expanded/Contracted alleles
col_names <- c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")
Events_df$Check_and_1 <- apply(Events_df[col_names], 1, function(x) if (x[1] %in% x[c(3:4)] & x[1] %in% x[c(5:6)]) {NA} else {x[1]})
Events_df$Check_and_2 <- apply(Events_df[col_names], 1, function(x) if (x[2] %in% x[c(3:4)] & x[2] %in% x[c(5:6)]) {NA} else {x[2]})
Events_df$Check_or_1 <- apply(Events_df[col_names], 1, function(x) if (x[1] %in% x[c(3:4)] | x[1] %in% x[c(5:6)]) {NA} else {x[1]})
Events_df$Check_or_2 <- apply(Events_df[col_names], 1, function(x) if (x[2] %in% x[c(3:4)] | x[2] %in% x[c(5:6)]) {NA} else {x[2]})
#     Analsyse the single events
Single_Evenets_df <- single_events(Events_df[Events_df$check_sum != 0,])
#     Analsyse the dounle events
Double_Events_df <- double_events(Events_df[Events_df$check_sum == 0,])
#     Combine all Expansion/Contraction Events
All_Events_df <- rbind(Single_Evenets_df, Double_Events_df)
All_Events_df <- merge(All_Events_df, loci[c("Call_ID", "Chr", "Gene", "Region")], by="Call_ID")
All_Events_df$Event_Size <- abs(All_Events_df$Event_Size)
All_Reads_df <- rbind(Stable_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], Poly_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], Events_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")])
#     Save files
trios_save_data(tcltk::tk_choose.dir(caption = "Select Save directory"), "NGC_Trios", All_Events_df, Poly_df, Stable_df, mdata)








write.table(trios_df, "D:/Work_Storage/MSSNG/Trios_MSSNG/MSSNG_Trios_TRIOFORMAT.tab", row.names = T, col.names = F, sep = "\t", quote = F)

head(All_Events_df)

SameElements <- function(a, b) return(identical(sort(a), sort(b)))
test_df <- Poly_df
test_df <- test_df[test_df$Allele1_Units_1 != test_df$Allele2_Units_1 & test_df$Allele1_Units_2 != test_df$Allele2_Units_2 & test_df$Allele1_Units_3 != test_df$Allele2_Units_3,]
test_df$check_test <- apply(test_df[c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], 1, function(x){
  x <- as.numeric(x)
  if (identical(sort(x[1:2]), sort(x[3:4])) & identical(sort(x[1:2]), sort(x[5:6]))){1} else {0}  #if ( (x[1] %in% x[3:4] & x[2] %in% x[3:4]) & (x[1] %in% x[5:6] & x[2] %in% x[5:6]) ){1} else {0}
  } )
out_test <- test_df[test_df$check_test == 1,]
agg_length <- aggregate(out_test$Call_ID, by=list("Call_ID"=out_test$Call_ID), length)
agg_max <- aggregate(out_test[c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], by=list("Call_ID"=out_test$Call_ID), max)
agg_max$max <- apply(agg_max[c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], 1, function(x) max(x))
agg_max$min <- apply(agg_max[c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], 1, function(x) min(x))
agg_max$diff <- agg_max$max - agg_max$min

test_final <- merge(agg_length, agg_max[c("Call_ID", "max", "min", "diff")])
test_final <- merge(test_final, loci[c("Call_ID", "Region", "Gene")])
test_final_filtered <- test_final[test_final$Call_ID %in% c("7.57927914", "2.181738159", "9.67340260", "22.45636629", "20.36531485", "2.132121699", "12.105629940", "19.590008"),]

write.table(test_final_filtered, "/home/dannear/Dropbox/University/PhD/Students/Kim/targets.xlsx", sep="\t", row.names = F)

#check records
col_inc <- c("Call_ID", "Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")
compound_df <- rbind(All_Events_df[col_inc], Poly_df[col_inc], Stable_df[col_inc])
nrow(unique(compound_df))
install.packages("sqldf")
library(sqldf)
a1NotIna2 <- sqldf('SELECT * FROM compound_df EXCEPT SELECT * FROM trios_df')
unique()
