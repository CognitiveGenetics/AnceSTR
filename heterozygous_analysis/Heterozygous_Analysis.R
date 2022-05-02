all_data <- informative_analysis(all_inform)
ngc_data <- informative_analysis(ngc_i)
trios_data <- informative_analysis(trio_i)
quads_data <- informative_analysis(quad_i)

mvf <- aggregate(list("count"=all_data$Transfer_Type), by=list("parent"=all_data$parent, "type"=all_data$Transfer_Type), length)
d1 <- ggplot(mvf, aes(x=parent, y=count)) + geom_bar(stat='identity', aes(fill=type), position = position_dodge(), width=.5)
d1 <- d1 + theme(legend.position="top", axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
d1 <- d1 + xlab("") + ylab("No. of Transfers\n") + scale_fill_discrete(name = "Allele Type", labels = c("Larger", "Smaller"))
d1

new_data <- trios_data
bias_diff <- aggregate(list("count"=new_data$Allele_Difference), by=list("type"=new_data$Transfer_Type, "diff"=new_data$Allele_Difference), length)
bias_diff <- merge(bias_diff[bias_diff$type == "S",], bias_diff[bias_diff$type == "L",], by="diff", all=T)
bias_diff$type.x <- "S"
bias_diff$type.y <- "L"
bias_diff[is.na(bias_diff)] <- 0
bias_diff$bias <- bias_diff$count.x/(bias_diff$count.x+bias_diff$count.y)*100
bias_diff$bias_t <- bias_diff$bias-50
bias_diff$ratio <- bias_diff$count.x/bias_diff$count.y

d2 <- ggplot(bias_diff, aes(x=diff, y=bias, label=bias)) + geom_point(stat = "identity", shape =2) + geom_line(stat = "identity")
d2 <- d2 + scale_x_continuous(name ="\n Allele Size Difference")#, breaks=seq(0,60, 10), limits = c(0, 60))
d2 <- d2 + scale_y_continuous(name ="Proportion Smaller Allele Selected (%)\n", breaks = seq(50, 100, 10), limits = c(50, 100))
d2 <- d2 + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
d2

d3 <- ggplot(bias_diff, aes(x=diff, y=ratio, label=ratio)) + geom_bar(stat = "identity")
d3 <- d3 + scale_x_continuous(name ="\n Allele Size Difference", breaks=seq(0,100, 10), limits = c(0, 100))
d3 <- d3 + scale_y_continuous(name ="Selection Ratio (Smaller Allele : Larger Allele) \n", breaks = seq(0, 20, 1), limits = c(0, 20))
d3 <- d3 + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
d3


build_d1 = function(all, ngc, trios, quads){
  for (x in c(myfunc(all), myfunc(ngc), myfunc(trios), myfunc(quads))){
    data <- get(x)
    assign(paste0("mvf_", x), aggregate(list("count"=data$Transfer_Type), by=list("parent"=data$parent, "type"=data$Transfer_Type), length))
    assign(paste0("d1_", x), ggplot(get(paste0("mvf_", x)), aes(x=parent, y=count)) + geom_bar(stat='identity', aes(fill=type), position = position_dodge(), width=.5) +
      theme(legend.position="top", axis.text.x = element_text(size = 13, face = "bold", colour = "black"), axis.text.y = element_text(size = 14, face = "plain", colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(), legend.text = element_text(colour="black", size = 14, face = "bold")) +
      xlab("") + ylab("No. of Transfers\n") + scale_fill_manual(name = "Allele Type", labels = c("Larger", "Smaller"), values = c("L" = "#0D0887FF", "S" = "#A62098FF"))
    )
  }
  grid.arrange(d1_all, d1_ngc, d1_trios, d1_quads)
}
build_d1(all_data, ngc_data, trios_data, quads_data)

build_d2 = function(all, ngc, trios, quads){
  for (x in c(myfunc(all), myfunc(ngc), myfunc(trios), myfunc(quads))){
    data <- get(x)
    bias_diff <- aggregate(list("count"=data$Allele_Difference), by=list("type"=data$Transfer_Type, "diff"=data$Allele_Difference), length)
    bias_diff <- merge(bias_diff[bias_diff$type == "S",], bias_diff[bias_diff$type == "L",], by="diff", all=T)
    bias_diff$type.x <- "S"
    bias_diff$type.y <- "L"
    bias_diff[is.na(bias_diff)] <- 0
    bias_diff$bias <- bias_diff$count.x/(bias_diff$count.x+bias_diff$count.y)*100
    bias_diff$bias_t <- bias_diff$bias-50
    bias_diff$ratio <- bias_diff$count.x/bias_diff$count.y
  
    assign(paste0("d2_", x), ggplot(bias_diff, aes(x=diff, y=bias, label=bias)) + geom_point(stat = "identity", shape=23, fill="black", color="red") + geom_text(label="Mendelian\nExpectation", x=15, y=45 , color = "Red") + #geom_line(stat = "identity") +
           scale_x_continuous(name ="\n Allele Size Difference", breaks=seq(0,100, 10), limits = c(0, 100)) +
           scale_y_continuous(name ="Proportion Smaller Allele Selected (%)\n", breaks = seq(40, 100, 10), limits = c(40, 100)) +
           theme(axis.line = element_line(size = 0.8, colour = "black", linetype=1), legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "black"), axis.text.y = element_text(size = 14, face = "bold", colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.background = element_blank(), legend.text = element_text(colour="black", size = 14, face = "bold")) +
           geom_hline(yintercept=50, linetype="dashed", color = "red") 
    )
  }
  grid.arrange(d2_all, d2_ngc, d2_trios, d2_quads)
}
build_d2(all_data, ngc_data, trios_data, quads_data)

build_d3 = function(all, ngc, trios, quads){
  for (x in c(myfunc(all), myfunc(ngc), myfunc(trios), myfunc(quads))){
    data <- get(x)
    bias_diff <- aggregate(list("count"=data$Allele_Difference), by=list("type"=data$Transfer_Type, "diff"=data$Allele_Difference), length)
    bias_diff <- merge(bias_diff[bias_diff$type == "S",], bias_diff[bias_diff$type == "L",], by="diff", all=T)
    bias_diff$type.x <- "S"
    bias_diff$type.y <- "L"
    bias_diff[is.na(bias_diff)] <- 0
    bias_diff$bias <- bias_diff$count.x/(bias_diff$count.x+bias_diff$count.y)*100
    bias_diff$bias_t <- bias_diff$bias-50
    bias_diff$ratio <- bias_diff$count.x/bias_diff$count.y
    bias_diff <- bias_diff[bias_diff$ratio != Inf,]
    
    assign(paste0("d3_", x),ggplot(bias_diff, aes(x=diff, y=ratio, label=ratio)) + geom_bar(stat = "identity", fill="Dark blue", width=1) + geom_text(label="Mendelian\nExpectation", x=95, y=2.5 ,label.size = 0.5, color = "Red") +
           scale_x_continuous(name ="\n Allele Size Difference", breaks=seq(0,100, 10), limits = c(0, 100)) +
           scale_y_continuous(name ="Selection Ratio (Smaller Allele : Larger Allele) \n", breaks = seq(0, 18, 2), limits = c(0, 18)) +
           theme(axis.line = element_line(size = 0.8, colour = "black", linetype=1), legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "black"), axis.text.y = element_text(size = 14, face = "plain", colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), panel.background = element_blank()) +
           geom_hline(yintercept = 1, colour="red", linetype = "longdash")
    )
  }
  grid.arrange(d3_all, d3_ngc, d3_trios, d3_quads)
}
build_d3(all_data, ngc_data, trios_data, quads_data)


bias_diff <- aggregate(list("count"=new_data$Allele_Difference), by=list("type"=new_data$Transfer_Type, "diff"=new_data$Allele_Difference, "Parent"=new_data$parent), length)
bias_diff <- merge(bias_diff[bias_diff$type == "S",], bias_diff[bias_diff$type == "L",], by=c("diff", "Parent"), all=T)
bias_diff <- bias_diff[c(1,2,4,6)]
colnames(bias_diff)[3] <- "smaller"
colnames(bias_diff)[4] <- "larger"
bias_diff$selection <- bias_diff$smaller/(bias_diff$smaller+bias_diff$larger)*100


d2 <- ggplot(bias_diff, aes(x=diff, y=selection, fill=Parent)) + geom_point(stat = "identity", shape =2) + geom_line(stat = "identity")
d2 <- d2 + scale_x_continuous(name ="\n Allele Size Difference")#, breaks=seq(0,60, 10), limits = c(0, 60))
d2 <- d2 + scale_y_continuous(name ="Proportion Smaller Allele Selected (%)\n", breaks = seq(50, 100, 10), limits = c(50, 100))
d2 <- d2 + theme(legend.position="top", legend.title=element_blank(), axis.text.x = element_text(size = 13, face = "bold", colour = "midnightblue"), axis.text.y = element_text(size = 14, face = "plain", colour = "midnightblue"), axis.title.x = element_text(size = 20, face = "bold", colour = "midnightblue"), axis.title.y = element_text(size = 20, face = "bold", colour = "midnightblue"), panel.background = element_blank(), legend.text = element_text(colour="midnightblue", size = 14, face = "bold"))
d2
