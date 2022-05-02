#     Function to get packages
Package_up = function(){
  list.of.packages <- c("stringr", "reshape2", "dplyr", "matrixStats", "zeallot", "ggplot2", "plyr", "ggsignif", "tidyr", "gridExtra")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if (length(new.packages)) install.packages(new.package)
  lapply(list.of.packages, require, character.only = TRUE)
}

#     Function to get mode
getmode = function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#     Function to determine largest/smallest expansion events (contraction or Expansion)
Get_Diff = function(x, minmax){
  if (minmax == "min"){
    largest <- 10000
    for (n in 1:length(x)){
      if (x[n] == 0) {next}
      if (abs(x[n]) <= abs(largest)){largest <- x[n]}
    } }
  if (minmax == "max"){
    largest <- 0
    for (n in 1:length(x)){
      if (abs(x[n]) >= abs(largest)){largest <- x[n]}
    } }  
  return(largest)
}

#     Function to calculate the most likely combination of repeat inheritence when two events occour at the same locus
calculate = function(sr){
  #     Format input data
  sr <- as.data.frame(t(sr))
  for (x in (1:length(sr[1,]))){
    if (x > 2){sr[,x] <- as.numeric(as.character(sr[,x]))}}
  #     Test the different combinations of parents and proband alleles
  trio <- sr[3:8]
  names(trio) <- c("p1", "p2", "m1", "m2", "f1", "f2")
  y <- data.frame(mother=c(trio$m1, trio$m1, trio$m2, trio$m2, trio$m1, trio$m1, trio$m2, trio$m2),
                  father=c(trio$f1, trio$f2, trio$f1, trio$f2, trio$f1, trio$f2, trio$f1, trio$f2),
                  proband1=c(trio$p1, trio$p1, trio$p1, trio$p1, trio$p2, trio$p2, trio$p2, trio$p2),
                  proband2=c(trio$p2, trio$p2, trio$p2, trio$p2, trio$p1, trio$p1, trio$p1, trio$p1))
  y$diff_m <- y$proband1 - y$mother
  y$diff_f <- y$proband2 - y$father
  #     Select for the most conservative change in repeat length
  y$sum <- apply(abs(y[c("diff_m", "diff_f")]), 1, sum)
  result <- y[y$sum==min(y$sum),]
  #     In the case that there is more than one outcome with the same expansion/contraction sizes, prioratize expansions.
  result$sum <- apply(result[c("diff_m", "diff_f")], 1, sum)
  result <- result[result$sum==max(result$sum),]
  #     Remove duplicated outputs
  result <- unique(result)
  #     If still more than 1 outcome parents share an allele of same size and which parent cannot be deterined. Else if only 1 out parent of orgin can be called.
  if (nrow(result) > 1){parent_m <- parent_f <- "Undetermined"} else {
    parent_m <- "Mother"
    parent_f <- "Father"
  }
  #     Detemine if result is a expansion or contration event
  if (result$diff_m[1] > 0){event_type_m <- "Expansion"} else {event_type_m <- "Contraction"}
  if (result$diff_f[1] > 0){event_type_f <- "Expansion"} else {event_type_f <- "Contraction"}  
  #     Return results
  return(list(result$proband1[1], result$diff_m[1], event_type_m, parent_m, result$proband2[1], result$diff_f[1], event_type_f, parent_f))
}

#     Function to import data
Import_Data = function(OS, build){
  if (build %in% c(19, 37)){if (OS == "Windows"){
    ngc <- read.csv(file="D:/Work_Storage/Cambridge_Data/NGC_37/CGG_Repeats_NGC_37_filtered.csv", header=TRUE, sep="\t")
    loci <- read.csv("D:/Work_Storage/Cambridge_Data/NGC_37/CGG_Repeats_NGC_37_filtered_Locus_summary.csv", header = TRUE, sep="\t")
  } else if (OS == "Linux") {
    ngc <- read.csv("/home/dannear/clusterhome/Cambridge_Data/raw_data/NGC_37/CGG_Repeats_NGC_37_filtered.csv", header=TRUE, sep="\t")
    loci <- read.csv("/home/dannear/clusterhome/Cambridge_Data/raw_data/NGC_37/NGC_37_filtered_Locus_summary.csv", header = TRUE, sep="\t")
  }} else
    if (build %in% c(38)){if (OS == "Windows"){
      ngc <- read.csv(file="D:/Work_Storage/Cambridge_Data/NGC_38/CGG_Repeats_NGC_38_filtered.csv", header=TRUE, sep="\t")
      loci <- read.csv("D:/Work_Storage/Cambridge_Data/NGC_38/NGC_38_Locus_summary.csv", header = TRUE, sep="\t")
    } else if (OS == "Linux") {
      ngc <- read.csv("/home/dannear/clusterhome/Cambridge_Data/Summaries/NGC_38/CGG_Repeats_NGC_38_filtered.csv", header=TRUE, sep="\t")
      loci <- read.csv("/home/dannear/clusterhome/Cambridge_Data/Summaries/NGC_38/NGC_38_Locus_summary.csv", header = TRUE, sep="\t")
    }}
  return(list(ngc, loci))  
}

#     Function to organise data into trios
Get_Trios = function(raw_df){
  #     IDENTIFY PATIENT TYPE  (01 = Proband, 02 = Mother, 03 = Father, 04+ = Sibling)
  if ("Member" %in% names(raw_df)) raw_df <- raw_df[c("Call_ID", "Member", "Sample_ID", "Allele1_Units", "Allele2_Units")] else {
    raw_df$Member <- as.numeric(str_split_fixed(raw_df$Sample_ID, "_0", 2)[,2])
    raw_df$Sample_ID <- str_split_fixed(raw_df$Sample_ID, "_0", 2)[,1]
    raw_df <- raw_df[c("Call_ID", "Member", "Sample_ID", "Allele1_Units", "Allele2_Units")]
  }
  #     REMOVE INCOMPLETE CASES
  auto_df <- raw_df[str_split_fixed(raw_df$Call_ID, "\\.", 2)[,1] != "X",]
  auto_df[auto_df == 0] <- NA
  sex_df <- raw_df[str_split_fixed(raw_df$Call_ID, "\\.", 2)[,1] == "X",]
  sex_df$Allele1_Units <- apply(sex_df, 1, function(x){ if (as.numeric(x[4]) == 0) NA else x[4]})
  raw_df <- rbind(auto_df, sex_df)
  raw_df <- raw_df[complete.cases(raw_df), ]
  #     GROUP DATA BY FAMILIES  
  fam <- melt(data = raw_df, id.vars = c("Call_ID", "Sample_ID", "Member"))
  fam$variable <- str_c(fam$variable, fam$Member, sep = "_")
  if (4 %in% unique(fam$member)){type <- 4} else {type <- 3}
  fam <- select(fam, -Member)
  fam <- dcast(fam, Call_ID + Sample_ID ~ variable, value.var = "value", drop=TRUE)
  if (type == 4) fam <- fam[,c("Call_ID", "Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3", "Allele1_Units_4", "Allele2_Units_4")]
  if (type == 3) fam <- fam[,c("Call_ID", "Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")]
  fam <- fam[complete.cases(fam), ]
  #     FILTER OUT SINGLES 
  fam_filtered <- fam[!(rowSums(is.na(fam[,3:(type*2)+2])) >= type*2),]
  #     FILTER OUT READS MISSING PROBANDS
  fam_filtered <- fam_filtered[!(rowSums(is.na(fam_filtered[,3:4])) == 2),]
  #     SELECT ALL TRIOS (Remove Siblings) and DUOs (Proband and parent)
  duos_trios <- fam_filtered[,c(1:8)]
  trios <- duos_trios[rowSums(is.na(duos_trios[,c(3:8)])) == 0,]

  return(list(trios, fam))
}

#     Function to collect data on expansion events and split data into stable, poly, and event loci
#     Set sizes (Trios, Quartets, Quintets)
Collect_Expansion_Events = function(df, set_size){
  if (set_size == "trios"){
    proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    ps1s2 <- proband} else
  if (set_size == "quartets"){
    proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    sibling1 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_4",  "Allele2_Units_4",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    sibling1$Sample_ID <- paste0(sibling1$Sample_ID, "_2")
    names(sibling1) <- names(proband)
    ps1s2 <- rbind(proband, sibling1)
  } else
  if (set_size == "quintets"){
    proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    sibling1 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_4",  "Allele2_Units_4",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    sibling2 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_5",  "Allele2_Units_5",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    sibling1$Sample_ID <- paste0(sibling1$Sample_ID, "_2")
    sibling2$Sample_ID <- paste0(sibling1$Sample_ID, "_3")
    names(sibling1) <- names(proband)
    names(sibling2) <- names(proband)
    ps1s2 <- rbind(proband, sibling1, sibling2)
  }

  
  ps1s2 <- ps1s2[complete.cases(ps1s2),]
  
  ps1s2$check_1_mother <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele1_Units_2, 1, 0)
  ps1s2$check_2_mother <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele2_Units_2, 1, 0)
  ps1s2$check_3_mother <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele1_Units_2, 1, 0)
  ps1s2$check_4_mother <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele2_Units_2, 1, 0)
  ps1s2$check_1_father <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele1_Units_3, 1, 0)
  ps1s2$check_2_father <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele2_Units_3, 1, 0)
  ps1s2$check_3_father <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele1_Units_3, 1, 0)
  ps1s2$check_4_father <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele2_Units_3, 1, 0)
  
  ps1s2$check_sum <- rowSums(ps1s2[,9:16]) 
  ps1s2$check_sum_1 <- rowSums(ps1s2[,c(9, 10, 13, 14)])
  ps1s2$check_sum_2 <- rowSums(ps1s2[,c(11, 12, 15, 16)])
  ps1s2$check_sum_mother <- rowSums(ps1s2[,c(9, 10, 11, 12)])
  ps1s2$check_sum_father <- rowSums(ps1s2[,c(13, 14, 15, 16)])
  
  stable_reps <- ps1s2[ps1s2$check_sum == 8, c(1:8, 17)]
  poly_reps <- ps1s2[ps1s2$check_sum %in% c(0, 1, 2, 3, 4, 6), c(1:8, 17)]
  event_reps <- ps1s2[(ps1s2$check_sum %in% c(0, 1)) | 
                        (ps1s2$check_sum == 3 & ps1s2$check_sum_1 %in% c(0,3) & ps1s2$check_sum_2 %in% c(0,3)) | 
                          (ps1s2$check_sum == 4 & ps1s2$check_sum_1 %in% c(0,4) & ps1s2$check_sum_2 %in% c(0,4)) | (ps1s2$check_sum == 4 & ps1s2$check_sum_mother %in% c(0,4) & ps1s2$check_sum_father %in% c(0,4)) | 
                            (ps1s2$check_sum == 2 & ps1s2$check_sum_1 == 1 & ps1s2$check_sum_2 == 1 & ps1s2$Allele1_Units_1 == ps1s2$Allele2_Units_1) | (ps1s2$check_sum == 2 & ps1s2$check_sum_1 %in% c(0,2) & ps1s2$check_sum_2 %in% c(0,2)), c(1:8, 17)]
  
  poly_reps <- anti_join(poly_reps, event_reps)
  
  #meta_data <- data.frame("No_Trios"=length(unique(ps1s2$Sample_ID)), "Total_Reads"=length(ps1s2$Call_ID), "Stable_Reads"=length(stable_reps$Call_ID), "Poly_Reads"=length(poly_reps$Call_ID), "Event_Reads"=length(event_reps$Call_ID))
  meta_data <- data.frame("Subset"=c("No_Trios", "Total_Reads", "Stable_Reads", "Poly_Reads", "Event_Reads"), "Count"=c(length(unique(ps1s2$Sample_ID)), length(ps1s2$Call_ID), length(stable_reps$Call_ID), length(poly_reps$Call_ID), length(event_reps$Call_ID)))
  
  return(list(event_reps, poly_reps, stable_reps, meta_data))
}
testCollect_Expansion_Events = function(df, set_size){
  if (set_size == "trios"){
    proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
    ps1s2 <- proband} else
      if (set_size == "quartets"){
        proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
        sibling1 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_4",  "Allele2_Units_4",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
        sibling1$Sample_ID <- paste0(sibling1$Sample_ID, "_2")
        names(sibling1) <- names(proband)
        ps1s2 <- rbind(proband, sibling1)
      } else
        if (set_size == "quintets"){
          proband <- df[c("Call_ID", "Sample_ID", "Allele1_Units_1",  "Allele2_Units_1",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
          sibling1 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_4",  "Allele2_Units_4",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
          sibling2 <- df[c("Call_ID", "Sample_ID", "Allele1_Units_5",  "Allele2_Units_5",  "Allele1_Units_2",  "Allele2_Units_2",  "Allele1_Units_3",  "Allele2_Units_3")]
          sibling1$Sample_ID <- paste0(sibling1$Sample_ID, "_2")
          sibling2$Sample_ID <- paste0(sibling1$Sample_ID, "_3")
          names(sibling1) <- names(proband)
          names(sibling2) <- names(proband)
          ps1s2 <- rbind(proband, sibling1, sibling2)
        }
  
  ps1s2 <- ps1s2[complete.cases(ps1s2),]
  
  ps1s2$check_1_mother <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele1_Units_2, 1, 0)
  ps1s2$check_2_mother <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele2_Units_2, 1, 0)
  ps1s2$check_3_mother <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele1_Units_2, 1, 0)
  ps1s2$check_4_mother <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele2_Units_2, 1, 0)
  ps1s2$check_1_father <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele1_Units_3, 1, 0)
  ps1s2$check_2_father <- ifelse(ps1s2$Allele1_Units_1 == ps1s2$Allele2_Units_3, 1, 0)
  ps1s2$check_3_father <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele1_Units_3, 1, 0)
  ps1s2$check_4_father <- ifelse(ps1s2$Allele2_Units_1 == ps1s2$Allele2_Units_3, 1, 0)
  
  ps1s2$check_X <- ifelse(as.character(str_split_fixed(ps1s2$Call_ID, "\\.", 2)[,1]) == "X", "X", "A")
  ps1s2$check_h <- ifelse(ps1s2$Allele1_Units_1 != ps1s2$Allele2_Units_1, 1, 0)
  ps1s2$check_sum <- rowSums(ps1s2[,9:16]) 
  ps1s2$check_sum_1 <- rowSums(ps1s2[,c(9, 10, 13, 14)])
  ps1s2$check_sum_2 <- rowSums(ps1s2[,c(11, 12, 15, 16)])
  ps1s2$check_sum_mother <- rowSums(ps1s2[,c(9, 10, 11, 12)])
  ps1s2$check_sum_father <- rowSums(ps1s2[,c(13, 14, 15, 16)])
  ps1s2$check_all <- paste0(ps1s2$check_X, ps1s2$check_h, ps1s2$check_sum_1, ps1s2$check_sum_2, ps1s2$check_sum_mother, ps1s2$check_sum_father)
  
  stable_reps <- ps1s2[ps1s2$check_all %in% c("A04444", "X13122",  "X03342"),]
  poly_reps <- ps1s2[ps1s2$check_all %in% c("A03342", "A03324", "A12112", "A11212", "A12121", "A11221", "A11111", "A13122", "A11322", "A12222", "A02222", "X12121", "X12112", "X11111", "X02222", "X11221"),]
  event_reps <- ps1s2[ps1s2$check_all %in% c("A10422", "A14022", "A13021", "A13012", "A10321", "A10312", "A12011", "A10101", "A12002", "A11010", "A11001", "A12020", "A10202", "A10220", "A10211", "A11120", "A10110", "A11102", "A10000", "A01120", "A01102", "A02240",  "A02204", "A00000", "X13021", "X10101", "X10321", "X12011", "X01120", "X11010", "X10000", "X11001", "X10110", "X10211", "X00000", "X12020", "X02240", "X11102", "X01102", "X10220"),]
  
  meta_data <- data.frame("Subset"=c("No_Trios", "Total_Reads", "Stable_Reads", "Poly_Reads", "Event_Reads"), "Count"=c(length(unique(ps1s2$Sample_ID)), length(ps1s2$Call_ID), length(stable_reps$Call_ID), length(poly_reps$Call_ID), length(event_reps$Call_ID)))
  
  return(list(event_reps, poly_reps, stable_reps, meta_data))
}

#     Function to analyse single expansion/contraction event loci
single_events = function(singles){
  #     Deterime which is the event allele
  singles$Event_Allele <- apply(singles, 1, function(x){ getmode(as.numeric(na.omit(x[25:28])))})
  singles <- singles[c(1:8, 29)]
  #     Compare event allele to parental alleles
  singles$Mother1 <- as.numeric(singles$Event_Allele) - as.numeric(singles$Allele1_Units_2)
  singles$Mother2 <- as.numeric(singles$Event_Allele) - as.numeric(singles$Allele2_Units_2)
  singles$Father1 <- as.numeric(singles$Event_Allele) - as.numeric(singles$Allele1_Units_3)
  singles$Father2 <- as.numeric(singles$Event_Allele) - as.numeric(singles$Allele2_Units_3)
  #     Determine the largest / smallest Expansion / Contraction Event
  singles$Event_Size_Max <- apply(singles[c("Mother1", "Mother2", "Father1", "Father2")], 1, function(x) Get_Diff(x, "max"))
  singles$Event_Size <- apply(singles[c("Mother1", "Mother2", "Father1", "Father2")], 1, function(x) Get_Diff(x, "min"))
  #     Classify Expansion or Contration
  singles$Event_Type <- sapply(singles$Event_Size, function(x) if (x < 0) {return("Contraction")} else if (x > 0) {return("Expansion")})
  singles$Parent <- ifelse((singles$Event_Size == singles$Mother1 | singles$Event_Size == singles$Mother2) & (singles$Event_Size == singles$Father1 | singles$Event_Size == singles$Father2), "Undetermined", 
                           ifelse((singles$Event_Size == singles$Mother1 | singles$Event_Size == singles$Mother2), "Mother", "Father"))
  singles <- singles[c("Call_ID", "Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3", "Event_Allele", "Event_Size", "Event_Type", "Parent")]
  #     Check for X chromosome for when proband is male
  x_singles <- filter(singles, Allele2_Units_1 == 0)
  singles <- filter(singles, Allele2_Units_1 != 0)
  x_singles$Parent <- "Mother"
  x_singles$Event_Size <- apply(x_singles[c("Allele1_Units_1", "Allele1_Units_2", "Allele2_Units_2")], 1, function(x) Get_Diff(c(as.numeric(x[1])-as.numeric(x[2]), as.numeric(x[1])-as.numeric(x[3])), "min"))
  singles <- rbind(singles, x_singles)
  return(singles)
}

#     Function to analyse double expansion/contraction event loci
double_events = function(doubles){
  doubles$result <- apply(doubles, 1, calculate)
  doubles <- suppressMessages(doubles %>% unnest_wider(result, names_sep=as.character(c(1:8)), names_repair = "minimal"))
  doubles_m <- doubles[c(1:8, 14:17)]
  doubles_f <- doubles[c(1:8, 18:21)]
  names(doubles_m) <- names(doubles_f) <- c("Call_ID", "Sample_ID", "Allele1_Units_1", "Allele2_Units_1", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3", "Event_Allele", "Event_Size", "Event_Type", "Parent")
  doubles_final <- rbind(doubles_m, doubles_f)
  doubles_final <- doubles_final[doubles_final$Event_Size != 0,]
return(doubles_final)
}

#     Function to deterime bias for smaller repeat sizes
heterozygous_parents = function(hetero, chr){
  # Define function to dfetermine if larger repeat allele is refelected in proband
  check_fun = function(x){
    if (x[3] != x[4] & x[5] == x[6]){
      if (max(x[3:4]) %in% x[1:2]) {"1_Mo"} else {"0_Mo"}
    } else
      if (x[3] == x[4] & x[5] != x[6]){
        if (max(x[5:6]) %in% x[1:2]) {"1_Fa"} else {"0_Fa"}
      } else
        if (x[3] != x[4] & x[5] != x[6]){  
          if (max(x[3:4]) %in% x[1:2] & min(x[5:6]) %in% x[1:2]) {"1_Mo-0_Fa"} else
            if (min(x[3:4]) %in% x[1:2] & max(x[5:6]) %in% x[1:2]) {"0_Mo-1_Fa"} else
              if (max(x[3:4]) %in% x[1:2] & max(x[5:6]) %in% x[1:2]) {"2_dbl"} else
                if (min(x[3:4]) %in% x[1:2] & min(x[5:6]) %in% x[1:2]) {"0_dbl"} else {NA}
        }
  }
  # Format data
  if (chr == "X") {hetero <- hetero[!( gsub("\\..*", "", hetero$Call_ID) == "X"),]}
  hetero <- hetero[hetero$Allele1_Units_2 != hetero$Allele2_Units_2 | hetero$Allele1_Units_3 != hetero$Allele2_Units_3, 1:8]
  col_names <- c("Call_ID", "Sample_ID", "proband1", "proband2", "parent1", "parent2")
  names(hetero) <- c("Call_ID", "Sample_ID", "proband1", "proband2", "mother1", "mother2", "father1", "father2")
  # Apply check function to data
  hetero$output<- apply(hetero[3:8], 1, function(x) check_fun(x))
  # Split data into destinct catagories, single events, double events where mother is recorded, and double events where father is recorded
  single <- hetero[complete.cases(hetero) & !(hetero$output %in% c("0_dbl", "2_dbl", "1_Mo-0_Fa", "0_Mo-1_Fa")),]
  double_Mo <- hetero[complete.cases(hetero) & hetero$output %in% c("0_dbl", "2_dbl", "1_Mo-0_Fa", "0_Mo-1_Fa"),]
  double_Fa <- hetero[complete.cases(hetero) & hetero$output %in% c("0_dbl", "2_dbl", "1_Mo-0_Fa", "0_Mo-1_Fa"),]
  double_Mo$output <- apply( double_Mo[9], 1, function(x) if (x %in% c("0_dbl", "0_Mo-1_Fa")) {"0_Mo"} else if (x %in% c("2_dbl", "1_Mo-0_Fa")){"1_Mo"} )
  double_Fa$output <- apply( double_Fa[9], 1, function(x) if (x %in% c("0_dbl", "1_Mo-0_Fa")) {"0_Fa"} else if (x %in% c("2_dbl", "0_Mo-1_Fa")){"1_Fa"} )
  # Bind organised events back together and spilt check and p_ID intop seperate columns
  hetero <- rbind(single, double_Mo, double_Fa)  
  hetero <- hetero %>% separate(output, c("check", "p_ID"), "_")
  # Deterime size of maximun repeat size present in each call
  hetero$max_mo <- apply(hetero[5:10], 1, function(x) if (x[6]=="Mo"){max(as.numeric(x[1:2]))} else  {NA} )
  hetero$max_fa <- apply(hetero[5:10], 1, function(x) if (x[6]=="Fa"){max(as.numeric(x[3:4]))} else  {NA} )
  # Melt data
  mhetero <- melt(hetero[c("p_ID", "check", "max_mo", "max_fa")], id=c("p_ID", "check"))
  mhetero <- mhetero[complete.cases(mhetero),]
  # Get the total number of times the larger allele was or was not transfered to child for both mother and father
  summ1 <- ddply(mhetero, .(p_ID, check), summarise, MAX = max(value), COUNT = length(value))
  # Get the total number of times the larger allele was or was not transfered to child for both mother and father by repeat size
  summ2 <- ddply(mhetero, .(p_ID, check, value), summarise, COUNT = length(value))
  # Get the total number of times the larger allele was or was not transfered to child ny repeat size
  summ3 <- ddply(mhetero, .(value, check), summarise, COUNT = length(value))
  
  return(list(summ1, summ2, summ3, mhetero))
}

#     Function to save raw data generated by trio analysis
trios_save_data = function(dir_path, save_name, All_Events_df, Poly_df, Stable_df, mdata){
  write.table(All_Events_df, paste0(dir_path, "/ancSTR_Trio_analysis_MutationReads_", save_name, ".tab"), sep="\t", row.names = F, col.names = T, quote=F)
  write.table(Poly_df, paste0(dir_path, "/ancSTR_Trio_analysis_InformativeReads_", save_name, ".tab"), sep="\t", row.names = F, col.names = T, quote=F)
  write.table(Stable_df, paste0(dir_path, "/ancSTR_Trio_analysis_MonogenicReads_", save_name, ".tab"), sep="\t", row.names = F, col.names = T, quote=F)
  write.table(mdata, paste0(dir_path, "/ancSTR_Trio_analysis_MetaData_", save_name, ".tab"), sep="\t", row.names = F, col.names = T, quote=F)
}

#     Function to format MSSNG input data (Proband = 1, Mother = 2, Father = 3, Sibling = 4)
assign_structure = function(reads_path, family_path, type){
  reads_df <- rename(read.csv(reads_path, sep="\t"), c("Sample_ID"="Index.ID"))
  family_manifest <- rename(read.csv(family_path, sep="\t"), c("Family.ID"="Sample_ID"))
  family_manifest <- family_manifest[family_manifest$Size==type,]
  samples <- unique(reads_df[c("Index.ID")])
  merge_df <- merge(x=samples, y=family_manifest, by="Index.ID")
  merge_df$Member <- apply(merge_df, 1, function(x){
    check <- toupper(as.character(x[3]))
    if (check %in% c("PROBAND", "AFFECTED_PROBAND", "PROBAND_AFFECTED")) 1 else
      if (check %in% c("MOTHER", "M", "MO")) 2 else
        if (check %in% c("FATHER", "F", "FA")) 3 else
           if (check %in% c("SIBLING", "UNAFFECTED_SIBLING", "AFFECTED_SIBLING", "SIBLING_UNAFFECTED", "SIBLING_AFFECTED")) 4})
  if (type == 4){
      check <- aggregate(list("sum"=merge_df$Member), by=list("Sample_ID"=merge_df$Sample_ID), sum)
      check <- merge(check[check$sum == 10,], merge_df, by="Sample_ID")
  } else 
    if (type == 3){
      check <- aggregate(list("sum"=merge_df$Member), by=list("Sample_ID"=merge_df$Sample_ID), sum)
      check <- merge(check[check$sum == 6,], merge_df, by="Sample_ID")
    }
  reads_df <- reads_df[reads_df$Index.ID %in% check$Index.ID,]
  reads_df <- merge(reads_df, merge_df[c("Index.ID", "Sample_ID", "Member")], by="Index.ID")
  return(unique(reads_df[c("Sample_ID", "Member", "Call_ID","Chr","Start","End","GT","Ref_Units","Allele1_Units","Allele2_Units")]))
}

#     Function to analyse informative reads
informative_analysis <- function(data){
  #filter out same mother == father genotypes
  data$s_check <- apply(data, 1, function(x) if (all(unique(x[c(5, 6)]) == unique(x[c(7, 8)]))){1} else {0} )
  data <- data[data$s_check != 1, c(1:8)]
  data_con <- data[data$s_check == 1, c(1:8)]
  #check which parents are heterozygous
  data$Chr <- str_split_fixed(data$Call_ID, "\\.", 2)[,1] # <- merge(data, loci[c("Call_ID", "Chr")], by="Call_ID")
  data$m_check <- apply(data[5:6], 1, function(x){ if (x[1] == x[2]){1} else 0})
  data$f_check <- apply(data[7:8], 1, function(x){ if (x[1] == x[2]){1} else 0})
  #sort event types
  data_b <- data[data$m_check == 0 & data$f_check == 0,]
  data_m <- data[data$m_check == 0 & data$f_check == 1,]
  data_m <- rbind(data_m, data_b)
  data_m$lists <- apply(data_m[c(3:8)], 1 , function(x){
    if (setequal(x[1:2], x[5:6])) {rem <- x[1:2]} else if (x[1] != x[2]){rem <- setdiff(x[1:2], x[5:6])} else {rem <- x[1:2]}
    if (min(x[3:4]) %in% rem){list(parent = "Mother", Allele_Size = min(x[3:4]), Transfer_Type = "S", Allele_Difference = abs(as.numeric(x[3])-as.numeric(x[4])))} else {list(parent = "Mother", Allele_Size = max(x[3:4]), Transfer_Type ="L", Allele_Difference = abs(as.numeric(x[3])-as.numeric(x[4])))} })
  data_m <- cbind(data_m[-c(9:12)], do.call(rbind.data.frame, data_m$lists))
  data_f <- data[data$m_check == 1 & data$f_check == 0,]
  data_f <- rbind(data_f, data_b)
  data_f <- data_f[data_f$Chr != "X",]
  data_f$lists <- apply(data_f[c(3:8)], 1 , function(x){
    if (setequal(x[1:2], x[3:4])) {rem <- x[1:2]} else if (x[1] != x[2]){rem <- setdiff(x[1:2], x[3:4])} else {rem <- x[1:2]}
    if (min(x[5:6]) %in% rem){list(parent = "Father", Allele_Size = min(x[5:6]), Transfer_Type = "S", Allele_Difference = abs(as.numeric(x[5])-as.numeric(x[6])))} else {list(parent = "Father", Allele_Size = max(x[5:6]), Transfer_Type ="L", Allele_Difference = abs(as.numeric(x[5])-as.numeric(x[6])))} })
  data_f <- cbind(data_f[-c(9:12)], do.call(rbind.data.frame, data_f$lists))
  data <- rbind(data_m, data_f)
  return(data)
}

#     Analysis functions
hetero_an <- function(Poly_df){
  c(data1, data2, data3, print_df) %<-% heterozygous_parents(Poly_df, "X")
  data3$Repeat_Bias <- apply(data3, 1, function(x) as.numeric(x[3])/sum(data3$COUNT[data3$value == as.numeric(x[1])])*100)
  data3$p_ID <- "All"
  data2$Repeat_Bias <- apply(data2, 1, function(x) as.numeric(x[4])/sum(data2$COUNT[data2$value == as.numeric(x[3]) & data2$p_ID == x[1]])*100)
  plot1 <- rbind(data2, data3)
  plot1$tot_Count <- apply(plot1, 1, function(x) sum(plot1$COUNT[plot1$value == as.numeric(x[3]) & plot1$p_ID == x[1]]))
  plot1 <- plot1[plot1$check==1,]
  
  b1 <- ggplot(data1, aes(x=p_ID, y=COUNT, fill=check)) + geom_bar(stat="identity", position=position_dodge())
  b1 <- b1 + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
  b1 <- b1 + scale_x_discrete(name ="\n Parent") 
  b1 <- b1 + scale_y_continuous(name="") + scale_fill_manual(values = c("deepskyblue", "red2"), labels=c("Smaller Allele", "Larger Allele"))
  
  
  #b2_plot <- aggregate(list("Total"=data2$COUNT), by=list("Larger.Repeat"=data2$value), FUN=sum)
  #b2 <- ggplot(b2_plot[b2_plot$Larger.Repeat > 19,], aes(x=Larger.Repeat, y=Total)) + geom_point() + geom_line() + theme_classic()
  #b2 <- b2 + scale_y_continuous(name="Total Count\n", breaks=seq(0,7000,1000), limits=c(0,7000)) 
  #print(b2)
  
  b3 <- ggplot(plot1[plot1$value>3,], aes(x=value, y=Repeat_Bias, group=p_ID)) + geom_line(aes(color=p_ID), size=1)  + geom_point(aes(color=p_ID)) + theme_classic() 
  b3 <- b3 + scale_x_continuous(name ="\n Larger Allele Repeat Size", breaks=seq(0,130, 10)) 
  b3 <- b3 + scale_y_continuous(name="% Selection Bias\n", breaks=seq(0,100,10), limits=c(0,100)) 
  b3 <- b3 + scale_color_manual(values=c('black', 'blue', 'deeppink'))
  b3 <- b3 + geom_hline(yintercept=50, linetype="dashed", color = "red") 
  b3 <- b3 + geom_vline(xintercept=50, linetype="dashed", color = "blue")
  b3 <- b3 + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
  b3 <- b3 + scale_colour_manual(labels=c("Total", "Father", "Mother"), values=c("midnightblue", "deepskyblue", "red"))
  print(b3)
  #apply(data2, 1, function(x) as.numeric(x[4])/sum(data2$COUNT[data2$value == as.numeric(x[3])])*100)
}
gene_an = function(All_Events_df, loci){
  genes <- All_Events_df#merge(All_Events_df, loci[c("Call_ID", "Gene", "Region")], by="Call_ID")
  mdata <- melt(genes[c("Call_ID", "Event_Type", "Gene", "Region", "Event_Size")], id=c("Call_ID", "Event_Type", "Gene", "Region"))
  agg <- ddply(mdata, .(Call_ID, Gene, Region, Event_Type), summarise,
               AVE = mean(value),
               SD = sd(value),
               VAR = var(value),
               MIN = min(value),
               MAX = max(value),
               LEN = length(value))
  agg$RANGE <- apply(agg[c("MIN","MAX")], 1, function(x) diff(x))
  g1 <- ggplot(agg[agg$Call_ID != "12.7781294",], aes(x=LEN, y=AVE)) + geom_point()
  print(g1)
  g2 <- ggplot(agg[agg$Event_Type == "Expansion" & agg$Call_ID != "12.7781294",], aes(x=LEN, y=RANGE)) + geom_point()
  print(g2)
  g3 <- ggplot(agg[agg$Event_Type == "Contraction" & agg$Call_ID != "12.7781294",], aes(x=LEN, y=RANGE)) + geom_point()
  print(g3)
  g4 <- ggplot(agg[agg$Call_ID != "12.7781294",], aes(x=LEN, y=RANGE)) + geom_point()
  print(g4)
  return(list(genes, mdata))
}
region_an = function(mdata, loci, genes){
  mdata <- mdata[abs(mdata$value) != 1,]
  genes$Event_Size <- apply(genes[c("Event_Allele", "Event_Size", "Event_Type")], 1, function(x) if (x[3] == "Contraction"){as.numeric(x[2])*-1} else if (x[3] == "Expansion"){as.numeric(x[2])} )
  region <- merge( ddply(mdata, .(Region, Event_Type), summarise, amount= length(value) ), aggregate(list("Total_length"=loci$Region), by=list("Region"=loci$Region), FUN=length), by="Region")
  region$Normalised <- region$amount/region$Total_length
  r1 <- ggplot(data=region, aes(x=Event_Type, y=Normalised, fill=Region)) + geom_bar(stat="identity", position="dodge")
  r1 <- r1 + scale_y_continuous(name="Ave. Devations\nper Repeat Locus\n", breaks=seq(0,120,20), limits=c(0,120)) + scale_x_discrete(name="")
  r1 <- r1 + theme(legend.position="top", axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.title = element_text(colour="midnightblue", size=16, face="bold"), legend.text = element_text(colour="black", size = 13, face = "bold"))
  r1 <- r1 + scale_fill_discrete(labels=c("Downstream", "Exonic", "Intergenic", "Intronic", "ncRNA", "Upstream", "3'-UTR", "5'-UTR"))
  r3 <- ggplot(data=genes[abs(genes$Event_Size) != 1,], aes(x=Region, y=Event_Size, fill=Region)) + geom_boxplot(outlier.shape = NA)
  r3 <- r3 + scale_y_continuous(name ="Mutation Size\n\n", breaks=seq(-30,40, 10), limits=c(-30,40))
  r3 <- r3 + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + geom_hline(yintercept=0, linetype="dashed", color = "Black") 
  r3 <- r3 + theme(legend.position='none', axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_blank(), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.title = element_text(colour="midnightblue", size=16, face="bold"), legend.text = element_text(colour="black", size = 13, face = "bold"))
  r3 <- r3 + scale_fill_discrete(labels=c("Downstream", "Exonic", "Intergenic", "Intronic", "ncRNA", "Upstream", "3'-UTR", "5'-UTR"))
  print(grid.arrange(r1, r3))
}
event_an = function(All_Events_df, genes){
  All_Events_df$Event_Size <- apply(All_Events_df, 1, function(x) if (x[11] == "Expansion") as.numeric(x[10]) else if (x[11] == "Contraction") as.numeric(x[10])*-1 )
  agg1 <- aggregate(All_Events_df$Event_Type, by = list(All_Events_df$Event_Type), length)
  agg2 <- aggregate(All_Events_df$Event_Size, by = list(All_Events_df$Event_Size), length)
  e1 <- ggplot(data=agg1, aes(x=Group.1, y=x)) + geom_bar(stat="identity") + geom_text(aes(label = x))
  e2 <- ggplot(data=agg2, aes(x=Group.1, y=x)) + geom_bar(stat="identity", colour="Dark red", fill="Red") 
  e2 <- e2  + theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.position="top", axis.line = element_line(size=1, colour="Black"), axis.ticks = element_line(size = 1), legend.title = element_text(face = "bold"),  panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() ) +
    #scale_y_continuous(name="Number of Mutations\n") + scale_x_continuous(name="\nMutation Size") #, breaks=seq(-20,60,20), limits=c(-20,60)

  write.table(agg2, "C:/Users/DLE/Dropbox/University/PhD/Bioinfo/agg2.tab", quote=F, col.names = T, row.names = F, sep = "\t")
  genes$X <- apply(genes[c("Event_Allele", "Event_Size", "Event_Type")], 1, function(x) if (x[3] == "Contraction"){as.numeric(x[1])+as.numeric(x[2])} else if (x[3] == "Expansion"){as.numeric(x[1])-as.numeric(x[2])} )   
  genes$Event_Size <- apply(genes[c("Event_Allele", "Event_Size", "Event_Type")], 1, function(x) if (x[3] == "Contraction"){as.numeric(x[2])*-1} else if (x[3] == "Expansion"){as.numeric(x[2])} )
  e3 <- ggplot(data=genes[genes$Call_ID != 	"12.7781294" & genes$X > 0,], aes(x=X, y=Event_Size, group=as.numeric(X))) + geom_boxplot(outlier.size = 0.3)
  e3 <- e3 + scale_x_continuous(name ="Parental Repeat Length\n " , breaks=seq(0,170, 10), limits=c(0, 170)) 
  e3 <- e3 + scale_y_continuous(name ="\nDeviation Size", breaks=seq(-100,250, 50),  limits=c(-100,200))
  e3 <- e3 + theme(axis.text.x = element_text(size = 13, face = "bold", colour = "Black"), axis.text.y = element_text(size = 14, face = "plain", colour = "Black"), axis.title.x = element_text(size = 20, face = "bold", colour = "Black"), axis.title.y = element_text(size = 20, face = "bold", colour = "Black"), panel.background = element_blank()) + coord_flip()
  print(e1)
  print(e2)
  print(e3)
}
parent_an =function(All_Events_df){
  p_agg <- aggregate(All_Events_df$Parent, by = list(All_Events_df$Parent), length)
  p1 <- ggplot(p_agg, aes(x=Group.1, y=x, fill=Group.1))+geom_bar(width = 1, stat = "identity") + theme(legend.position="top")
  print(p1)
}

#     Function to get name of 
myfunc <- function(v1) {
  deparse(substitute(v1))
}


