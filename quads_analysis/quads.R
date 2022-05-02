("D:/Work_Storage/MSSNG/Quads_MSSNG/CGG_Repeats_Quads_Full_20220215_filtered.csv")
ngc <- assign_structure("D:/Work_Storage/MSSNG/Quads_MSSNG/CGG_Repeats_Quads_Full_20220215_filtered.csv", "C:/Users/DLE/Dropbox/University/PhD/Bioinfo/MSSNG/Sample Data/MSSNG_family_manifest.tab", 4)

probands <- ngc[ngc$Member != 4,]
sibling <- ngc[ngc$Member != 1,]
sibling[ngc$Member == 4] <- 1
sibling["Member"][sibling["Member"] == 4] <- 1

c(test_df, families_test) %<-% Get_Trios(probands)

Process_Trios = function(ngc){
c(trios, families) %<-% Get_Trios(ngc)
#     ORGANISE AND CLASSIFY TRIO READS
c(Events_df, Poly_df, Stable_df, mdata) %<-% testCollect_Expansion_Events(trios, "trios")
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
return(list(Events_df, Poly_df, Stable_df, mdata, All_Events_df, All_Reads_df))
}
c(pro_Events, pro_Info, pro_mono, pro_mdata, pro_All_Event, pro_All_Reads) %<-% Process_Trios(probands)
c(sib_Events, sib_Info, sib_mono, sib_mdata, sib_All_Event, sib_All_Reads) %<-% Process_Trios(sibling)
pro_All_Event$Type <- "Proband"
sib_All_Event$Type <- "Sibling"
quad_All_Events <- rbind(pro_All_Event, sib_All_Event)


event_counts <- aggregate(list("count"=quad_All_Events$Event_Size), by=list("Gene"=quad_All_Events$Gene, "Type"=quad_All_Events$Type), length)


test <- aov(count ~ Type * Gene, data = event_counts)
ttt <- TukeyHSD(test)

summary(quad_All_Events)

one.way <- aov(Event_Size ~ Type, data = quad_All_Events)
summary(one.way)

two.way <- aov(Event_Size ~ Gene + Type, data = quad_All_Events)
summary(two.way)

two.way2 <- aov(Event_Size ~ Gene * Type, data = quad_All_Events)
summary(two.way2)
tt <- TukeyHSD(two.way2)

length(pro_All_Event$Gene[pro_All_Event$Gene %in% genelist])
length(sib_All_Event$Gene[sib_All_Event$Gene %in% genelist])

test <- merge(aggregate(list("count"=pro_All_Event$Gene), by=list("gene"=pro_All_Event$Gene), length), aggregate(list("count"=sib_All_Event$Gene), by=list("gene"=sib_All_Event$Gene), length), by="gene")
test$diff <- test$count.x-test$count.y
test$ration <- test$count.x/test$count.y
test$ration2 <- apply(test, 1, function(x) if (as.numeric(x[2]) > as.numeric(x[3])) as.numeric(x[2])/as.numeric(x[3]) else if (as.numeric(x[2]) < as.numeric(x[3])) -1*as.numeric(x[3])/as.numeric(x[2]) else 0 )
test$length <- test$count.x+test$count.y
test$percentage.diff <- abs(test$diff)/test$length*100

test2 <- test[test$diff != 0,]

p<-ggplot(data=test, aes(x=gene, y=ration2)) +
  geom_bar(stat="identity")
p
p <- ggplot(quad_All_Events[quad_All_Events$Event_Size != 1,], aes(x=Type, y=Event_Size, fill=Gene)) + 
  geom_boxplot()
p <- p + theme(legend.position="none") 
p

lencheck <- merge(aggregate(list("total"=as.numeric(pro_All_Event$Event_Size)), by=list("gene"=pro_All_Event$Gene), sum), aggregate(list("total"=as.numeric(sib_All_Event$Event_Size)), by=list("gene"=sib_All_Event$Gene), sum), by="gene")





# Differences Probnad versus sibiling
com_melt1 <- melt(pro_All_Event[pro_All_Event$Event_Type == "Expansion", c("Gene", "Type", "Event_Size")], id=c("Gene", "Type"))
test1 <- ddply(com_melt1, .(Gene, Type), summarise, AVE = mean(value, na.rm=TRUE), LEN = length(value))
com_melt2 <- melt(sib_All_Event[sib_All_Event$Event_Type == "Expansion", c("Gene", "Type", "Event_Size")], id=c("Gene", "Type"))
test2 <- ddply(com_melt2, .(Gene, Type), summarise, AVE = mean(value, na.rm=TRUE), LEN = length(value))
test_merge <- merge(test1, test2, by = "Gene", all=T)
test_merge$Type.x <- 1 #proband
test_merge$Type.y <- 0 #sibling
test_bind <-  rbind(setNames(test_merge[c("Gene", "Type.x", "AVE.x", "LEN.x")], c("Gene", "Type", "AVE", "LEN")), setNames(test_merge[c("Gene", "Type.y", "AVE.y", "LEN.y")], c("Gene", "Type", "AVE", "LEN")))
test_bind[is.na(test_bind)] <- 0

model <- glm(Type ~ LEN + Gene, data = test_bind, family = binomial)
summary(model)
event_stats <- data.frame(summary(model)$coefficients)
names(event_stats) <- c("Estimate", "Std.Error", "z.value", "p.value")
event_stats$Gene <- str_replace(as.character(row.names(event_stats)), "Gene", "")
event_stats$sig.Size <- apply(event_stats, 1, function(x) if (as.numeric(x[4]) < 0.001) "***" else if (as.numeric(x[4]) < 0.01) "**" else if (as.numeric(x[4]) < 0.05) "*" else if (as.numeric(x[4]) <= 0.1) "." else "")

model_size <- glm(Type ~ LEN + AVE + Gene, data = test_bind, family = binomial)
size_stats <- data.frame(summary(model_size)$coefficients)
names(size_stats) <- c("Estimate", "Std.Error", "z.value", "p.value")
size_stats$Gene <- str_replace(as.character(row.names(size_stats)), "Gene", "")
size_stats$sig.Size <- apply(size_stats, 1, function(x) if (as.numeric(x[4]) < 0.001) "***" else if (as.numeric(x[4]) < 0.01) "**" else if (as.numeric(x[4]) < 0.05) "*" else if (as.numeric(x[4]) <= 0.1) "." else "")

