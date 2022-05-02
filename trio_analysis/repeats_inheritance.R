#     INTERGENERATIONAL REPEATS ANALYSIS
#     IMPORT FUNCTIONS
source("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/anceSTR/anceSTR/functions/functions.R")
#     LIBRARIES
Package_up()
#     IMPORT BASE DATA 1
c(ngc, loci) %<-% Import_Data(as.character(Sys.info()[['sysname']]), 38)
#     IMPORT BASE DATA 2
ngc <- read.csv(tcltk::tk_choose.files(caption = "Select Reads File"), header = T, sep = "\t")
loci <- read.csv(tcltk::tk_choose.files(caption = "Select Corresonding Loci File"), header = T, sep = "\t")
#     IMPORT BASE DATA 3
ngc <- assign_structure("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Analysis/Quads/CGG_Repeats_Quads_Full_20220215_filtered.csv","C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Sample Data/MSSNG_family_manifest.tab", 4)
#     ORGANISE RAW DATA INTO TRIOS
c(trios_df, families_df) %<-% Get_Trios(ngc)
#     ORGANISE AND CLASSIFY TRIO READS
c(Events_df, Poly_df, Stable_df, mdata) %<-% Collect_Expansion_Events(all_trios, "trios")
c(Events_df, Poly_df, Stable_df, mdata) %<-% testCollect_Expansion_Events(trio_t, "trios")
#     ANALYSE REPEAT LENGTH CHANGES
col_names <- c("Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")
Events_df$Check_and_1 <- apply(Events_df[col_names], 1, function(x) if (x[1] %in% x[c(3:4)] & x[1] %in% x[c(5:6)]) {NA} else {x[1]})
Events_df$Check_and_2 <- apply(Events_df[col_names], 1, function(x) if (x[2] %in% x[c(3:4)] & x[2] %in% x[c(5:6)]) {NA} else {x[2]})
Events_df$Check_or_1 <- apply(Events_df[col_names], 1, function(x) if (x[1] %in% x[c(3:4)] | x[1] %in% x[c(5:6)]) {NA} else {x[1]})
Events_df$Check_or_2 <- apply(Events_df[col_names], 1, function(x) if (x[2] %in% x[c(3:4)] | x[2] %in% x[c(5:6)]) {NA} else {x[2]})
#     Analsyse the single events
Single_Evenets_df <- single_events(Events_df[Events_df$check_all %in% c("A10422", "A14022", "A13021", "A13012", "A10321", "A10312", "A12011", "A10101", "A12002", "A11010", "A11001", "A12020", "A10202","A10220","A10211", "A11120","A10110", "A11102", "A01120", "A01102", "A02240", "A02204", "X13021", "X10101", "X10321", "X12011", "X01120", "X11010", "X10000", "X11001", "X10110", "X10211", "X02240", "X11102", "X01102", "X10220" ),])
#     Analsyse the dounle events
Double_Events_df <- double_events(Events_df[Events_df$check_all %in% c("X00000", "X12020", "A00000", "A10000"),c("Call_ID","Sample_ID","Allele1_Units_1","Allele2_Units_1","Allele1_Units_2","Allele2_Units_2","Allele1_Units_3","Allele2_Units_3","check_sum","Check_and_1","Check_and_2","Check_or_1","Check_or_2")])
#     Combine all Expansion/Contraction Events
All_Events_df <- rbind(Single_Evenets_df, Double_Events_df)
All_Events_df <- merge(All_Events_df, loci[c("Call_ID", "Chr", "Gene", "Region")], by="Call_ID")
All_Events_df$Event_Size <- abs(All_Events_df$Event_Size)
All_Reads_df <- rbind(Stable_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], Poly_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], Events_df[c("Call_ID","Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")])
#     Save files
trios_save_data(tcltk::tk_choose.dir(caption = "Select Save directory"), "MSSNG_Trios", All_Events_df, Poly_df, Stable_df, mdata)

"XYLT1", "FMR1","AFF2","HOXD1","LOC642361","LRP12","GIPC1", "NOTCH2NLC", "DIP2B", "AFF3", "FRA10AC1"
