
#     Import gene lists of dissease-associated genes
nddgenes <- read.csv("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/anceSTR/anceSTR/Data/Repeats_of_interest_NDD_ID.txt", header = T, sep="\t")
nddgeneslist <- unique(read.csv("C:/Users/DLE/Dropbox/University/PhD/Bioinfo/anceSTR/anceSTR/Data/Repeats_of_interest_NDD_ID_2.txt", header = T, sep="\t")[1])
fs_37_sites <- data.frame(FragileSite=c("FRA2A", "FRA7A", "FRA10A", "FRA11A", "FRA11B", "FRA12A", "FRA16A", "FRA19B", "FRAXA", "FRAXE", "FRAXF"), ID=c("2.100721250", "7.55955533", "10.95462280", "11.66512291", "11.119077000", "12.50898785", "16.17564765", "19.2308194", "X.146993568", "X.147582126", "X.148713408") )
fs_38_sites <- data.frame(FragileSite=c("FRA2A", "FRA7A", "FRA10A", "FRA11A", "FRA11B", "FRA12A", "FRA16A", "FRA19B", "FRAXA", "FRAXE", "FRAXF"), ID=c("2.100104473", "7.55887840", "10.93702523", "11.66744820", "11.119206290", "12.50505002", "16.17470908", "19.2308195", "X.147912050", "X.148500606", "X.149631714") )
#     Melt data present in trio format
trio_melt <- melt(all_trios[c("Call_ID", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], id="Call_ID")
#     Remove all 0 values, i.e. missing X chromosome alleles in males
trio_melt <- trio_melt[trio_melt$value != 0,]
#     Melt data of all event reads
event_melt <- melt(all_mutations[c("Call_ID", "Event_Size", "Event_Type", "Parent")], id=c("Call_ID", "Event_Type", "Parent"))
event_melt <- event_melt[event_melt$value != 0,]
#     Calculate event metrics per locus
ave_size_event <- ddply(trio_melt, .(Call_ID), summarise, AVE = mean(value), MED = median(value), Len = length(value)/2)
size_len_event <- ddply(event_melt, .(Call_ID), summarise, Size = mean(abs(value)), total = length(value))
type_len_event <- ddply(event_melt, .(Call_ID, Event_Type), summarise, total = length(value))
parent_len_event <- ddply(event_melt, .(Call_ID, Parent), summarise, total = length(value))
#     Merge data for plotting
tplot <- Reduce(function(x,y) merge(x = x, y = y, by = "Call_ID"), list(ave_size_event, size_len_event, loci[c("Call_ID", "Gene", "Region")]))
tplot$rate <- tplot$total/tplot$Len*100
tplot$rate2 <- tplot$total/tplot$Len
tplot$Gene_1 <- str_split_fixed(tplot$Gene, "->", 2)[,1]
tplot$Gene_2 <- str_split_fixed(tplot$Gene, "->", 2)[,2]
labplot <- tplot[tplot$Call_ID %in% fs_38_sites$ID,]
#     Build plot of mutation rate verses mean parental repeat length
aveplot <- ggplot(tplot, aes(x=AVE, y=rate2, fill=Region, color=Region, label=Gene)) + geom_point() +
  theme(legend.position="top", axis.line = element_line(size=1.5, colour="Black"), axis.ticks = element_line(size = 1), legend.title = element_text(face = "bold"), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE), name="Deviation Rate\n(Change per gamete per generation)\n", breaks=seq(0,1,0.1), limits=c(0,1)) + scale_x_continuous(name="\nMean Parental Repeat Length", breaks=seq(0,50,5), limits=c(0,50)) +
  theme(legend.title = element_text(size = 13, face = "bold", colour = "Black"), axis.text.x = element_text(size = 13, face = "bold", colour = "Black"), axis.text.y = element_text(size = 14, face = "plain", colour = "Black"), axis.title.x = element_text(size = 20, face = "bold", colour = "Black"), axis.title.y = element_text(size = 20, face = "bold", colour = "Black"), panel.background = element_blank(), legend.text = element_text(colour="Black", size = 13), axis.line = element_line(size = 1, colour = "Black", linetype=1))
aveplot
#   Select repeats that overlap with disease-associated genes
geneplot <- tplot[(tplot$Gene_1 %in% nddgenes$Gene |  tplot$Gene_2 %in% nddgenes$Gene)& tplot$AVE<60,]
aveplot2 <- aveplot + geom_point(data=geneplot, shape=25, color="red", size=3)+
  geom_smooth(method=lm, se=F,  aes(group=1)) +
  geom_text(data=subset(tplot, Call_ID %in% fs_38_sites$ID & Gene != "LINGO3" & Gene != "CBL"),  position=position_jitter(width=0.1,height=0.1), aes(label=Gene), hjust=-0.2, vjust=0, colour="Black", fontface=2) +
  geom_point(data=subset(tplot, Call_ID %in% fs_38_sites$ID), shape=21, color="Black", size=2)
aveplot2
#     Build plot of Ave mutation size verses mean parental repeat length
tplot$col1 <- "black" 
sizeplot <- ggplot(tplot, aes(x=AVE, y=Size, fill=Region, color=Region, label=Gene)) + geom_point() + stat_function(fun = function(x) {-0.3*x + 12}, size = 1, colour = "blue") +
  geom_point(data=subset(tplot, Call_ID %in% fs_38_sites$ID), shape=21, color="Black", size=2) +
  geom_text(data=subset(tplot, Call_ID %in% fs_38_sites$ID & Gene != "LINGO3" & Gene != "CBL" ),  position=position_jitter(width=2,height=2), aes(label=Gene), hjust=-0.2, vjust=0, colour="Black", fontface=2) +
  theme(legend.position="top", axis.line = element_line(size=1.5, colour="Black"), axis.ticks = element_line(size = 1), legend.title = element_text(face = "bold"), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(legend.title = element_text(size = 13, face = "bold", colour = "Black"), axis.text.x = element_text(size = 13, face = "bold", colour = "Black"), axis.text.y = element_text(size = 14, face = "plain", colour = "Black"), axis.title.x = element_text(size = 20, face = "bold", colour = "Black"), axis.title.y = element_text(size = 20, face = "bold", colour = "Black"), panel.background = element_blank(), legend.text = element_text(colour="Black", size = 13), axis.line = element_line(size = 1, colour = "Black", linetype=1)) +
  scale_y_continuous(name="Mean Deviation Size\n (CGG Units)\n", breaks=seq(0,40,5), limits=c(0,40)) + scale_x_continuous(name="\nMean Parental Repeat Length", breaks=seq(0,50,5), limits=c(0,50)) +
  scale_fill_discrete(labels=c("downstream"="Downstream", "exonic"="Exonic", "intergenic"="Intergenic", "intronic"="Intronic", "ncRNA"="ncRNA", "upstream"="Upstream", "UTR3"="3'-UTR", "UTR5"="5'-UTR")) 
sizeplot

#     Graphs for larger mutations outside the normal range
large <- tplot[tplot$rate2 > 0.01,]
rplot_rate <- ggplot(large, aes(x=Region, y=rate2, fill=Region)) + geom_boxplot() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE),name ="Deviation Rate\n(Events per gamete per generation)\n", breaks=seq(0, 1, 0.1), limits=c(0,1)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.position='top', axis.text.x = element_text(size = 13, face = "bold", colour = "Black"), axis.text.y = element_text(size = 14, face = "plain", colour = "Black"), axis.title.x = element_text(size = 20, face = "bold", colour = "Black"), axis.title.y = element_text(size = 20, face = "bold", colour = "Black"), panel.background = element_blank(), legend.title = element_text(colour="Black", size=16, face="bold"), legend.text = element_text(colour="Black", size = 13, face = "bold")) +
  scale_fill_discrete(labels=c("downstream"="Downstream", "exonic"="Exonic", "intergenic"="Intergenic", "intronic"="Intronic", "ncRNA"="ncRNA", "upstream"="Upstream", "UTR3"="3'-UTR", "UTR5"="5'-UTR")) +
  scale_x_discrete(name = "\nGenetic Region", labels=c("downstream"="Downstream", "exonic"="Exonic", "intergenic"="Intergenic", "intronic"="Intronic", "ncRNA"="ncRNA", "upstream"="Upstream", "UTR3"="3'-UTR", "UTR5"="5'-UTR"))
rplot_rate

rplot_size <- ggplot(large, aes(x=Region, y=Size, fill=Region)) + geom_boxplot() +
  scale_y_continuous(name ="Average Deviation Size\n(CGG Units)\n", breaks=seq(0, 35, 5), limits=c(0,35)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.position='top', axis.text.x = element_text(size = 13, face = "bold", colour = "Black"), axis.text.y = element_text(size = 14, face = "plain", colour = "Black"), axis.title.x = element_text(size = 20, face = "bold", colour = "Black"), axis.title.y = element_text(size = 20, face = "bold", colour = "Black"), panel.background = element_blank(), legend.title = element_text(colour="Black", size=16, face="bold"), legend.text = element_text(colour="Black", size = 13, face = "bold")) +
  scale_fill_discrete(labels=c("downstream"="Downstream", "exonic"="Exonic", "intergenic"="Intergenic", "intronic"="Intronic", "ncRNA"="ncRNA", "upstream"="Upstream", "UTR3"="3'-UTR", "UTR5"="5'-UTR")) +
  scale_x_discrete(name = "\nGenetic Region", labels=c("downstream"="Downstream", "exonic"="Exonic", "intergenic"="Intergenic", "intronic"="Intronic", "ncRNA"="ncRNA", "upstream"="Upstream", "UTR3"="3'-UTR", "UTR5"="5'-UTR"))
rplot_rate

sig_genes <- read.csv("C:/Users/DLE/Dropbox/University/PhD/Writing/Manuscript 2/Data analysis/scripts/repeat_stats.txt", sep="t")
siglist <- data.frame("genes"=unique(c(str_split_fixed(sig_genes$Gene, "->", 2)[,1], str_split_fixed(sig_genes$Gene, "->", 2)[,2])))
siglist <- siglist[siglist$genes %in% nddgeneslist$Gene,]


names(vari) <- c("x", "y")
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
p1 <- ggplot(vari, aes(x=x, y=y)) + geom_point() + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
p1 <- p1 + geom_text(label = lm_eqn(vari), parse = TRUE)
p1

vari2 <- tplot[c("rate2", "AVE")]
names(vari2) <- c("y", "x")
p2 <- ggplot(vari2, aes(x=x, y=y)) + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
p2 <- p2 + geom_text(label = lm_eqn(vari), parse = TRUE)
p2


