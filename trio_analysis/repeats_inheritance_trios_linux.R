#     INTERGENERATIONAL REPEATS ANALYSIS
#     IMPORT FUNCTIONS
source("/home/dannear/MSSNG_Data/scripts/functions.R")
#     LIBRARIES
Package_up()
#     IMPORT BASE DATA
loci <- read.csv("/home/dannear/MSSNG_Data/Trios/Trios_20220215_by_locus.csv", header = T, sep = "\t")
ngc <- assign_structure("/home/dannear/MSSNG_Data/Trios/CGG_Repeats_Trios_filtered.csv","/home/dannear/MSSNG_Data/Manifests/MSSNG_family_manifest.tab", 3)
#     ORGANISE RAW DATA INTO TRIOS
c(trios_df, families_df) %<-% Get_Trios(ngc)
#     ORGANISE AND CLASSIFY TRIO READS
c(Events_df, Poly_df, Stable_df, mdata) %<-% Collect_Expansion_Events(families_df, "quartets")
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
trios_save_data("/home/dannear/MSSNG_Data/Trios", "MSSNG_Trios", All_Events_df, Poly_df, Stable_df, mdata)
