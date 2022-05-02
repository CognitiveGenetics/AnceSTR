# Import data
library(circlize)

# Set data
plot_data <- melt(all_trios[c("Call_ID", "Allele1_Units_2", "Allele2_Units_2", "Allele1_Units_3", "Allele2_Units_3")], id="Call_ID")
plot_data <- plot_data[plot_data$value != 0,]
event_data <- melt(all_mutations[c("Call_ID", "Event_Size", "Event_Type")], id=c("Call_ID", "Event_Type"))
event_data <- event_data[event_data$value != 0,]

ave_size <- ddply(plot_data, .(Call_ID), summarise, AVE = mean(value), MED = median(value), Len = length(value)/2)
ave_size_2 <- ddply(event_data, .(Call_ID, Event_Type), summarise, AVE = mean(value), MED = median(value), Len = length(value))
ave_size_2$AVE <- apply(ave_size_2[c("Event_Type", "AVE")], 1, function(x) if (x[1] == "Contraction") -1*as.numeric(x[2]) else if (x[1] == "Expansion") as.numeric(x[2]) )
ave_size_2$start <- as.numeric(str_split_fixed(ave_size_2$Call_ID, "\\.", 2)[,2])
ave_size_2$col <- apply(ave_size_2[c("Event_Type", "AVE")], 1, function(x) if (x[1] == "Contraction") "#0D0887FF" else if (x[1] == "Expansion") "#A62098FF" )
ave_size_2 <- ave_size_2[ave_size_2$AVE > -20 & ave_size_2$AVE < 50,]
size_len <- ddply(event_data, .(Call_ID), summarise, Ave_Mut_Size = mean(abs(value)), total = length(value))
type_len <- ddply(event_data, .(Call_ID, Event_Type), summarise, event_count = length(value))

all_data <- Reduce(function(x, y) merge(x, y, all=TRUE, by="Call_ID"), list(ave_size, size_len, type_len))
all_data$Rate <- all_data$total/all_data$Len*100
all_data$chr <- paste0("chr", str_split_fixed(all_data$Call_ID, "\\.", 2)[,1])
all_data$start <- as.numeric(str_split_fixed(all_data$Call_ID, "\\.", 2)[,2])
mut_count <- melt(data = all_data[c("Call_ID", "Event_Type", "event_count")], id.vars = c("Call_ID", "Event_Type"))
mut_count <- select(mut_count, -variable)
mut_count <- dcast(mut_count, Call_ID ~ Event_Type, value.var = "value", drop=TRUE)[c("Call_ID", "Expansion", "Contraction")]
mut_count[is.na(mut_count)] <- 0
final <- merge(unique(all_data[c("Call_ID", "chr", "start", "AVE", "MED", "Len", "Ave_Mut_Size", "total", "Rate")]), mut_count, by="Call_ID")
final[is.na(final)] <- 0
final$colour_rate <- apply(final[c("Rate")], 1, function(x) if (x[1] > 30) "Red" else if (x[1] > 1) "Orange" else if  (x[1] > 0) "Green" else if (x[1] == 0) "Blue")
final$Contraction <- final$Contraction*-1
ave_size_2$chr <- paste0("chr", str_split_fixed(ave_size_2$Call_ID, "\\.", 2)[,1])
bed_labels <- merge(final[c("Call_ID", "chr", "start", "Rate")], loci[c("Call_ID", "Gene")])
bed_labels$end <- bed_labels$start+1
plot_labels <- bed_labels[bed_labels$Rate > 20,]
plot_labels <- plot_labels[c("chr", "start", "end", "Gene")]
plot_labels <- plot_labels[!duplicated(plot_labels$Gene), ]
clean_genes = function(pl){
  pl$Gene_1 <- str_split_fixed(pl$Gene, "->", 2)[,1]
  pl$Gene_2 <- str_split_fixed(pl$Gene, "->", 2)[,2]
  pl1 <- pl[c("chr","start","end","Gene")]
  pl1$Gene <- pl$Gene_1
  pl2 <- pl[c("chr","start","end","Gene")]
  pl2$Gene <- pl$Gene_2
  pl2$start <-  pl2$start + 100
  pl2$end <-  pl2$start + 1
  pl <- rbind(pl1, pl2)
  pl <- pl[pl$Gene != "" & pl$Gene != "NONE",]
  pl <- pl[!duplicated(pl$Gene), ]
  pl$Gene <- apply(pl[c("Gene")], 1, function(x) if (nchar(x[1])>10) substr(x[1],1,nchar(x[1])-3) else x[1])
  pl$cols <- apply(pl[c("chr")], 1, function(x){
    chr <- gsub("chr", "", x[1])
    if (chr == "X") "cyan4" else if  ((as.numeric(chr) %% 3) == 0) "cyan4" else if  ((as.numeric(chr) %% 2) == 0) "#0D0887FF" else "#A62098FF"} )
  return(pl)
}
plot_labels <- clean_genes(plot_labels)

#     PLOT LABELS
tiff("/home/dannear/Dropbox/University/PhD/Writing/Manuscript 2/Data analysis/graphs/circos/test2.tiff", units="in", width=8, height=8, res=300)
circos.clear()
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(plotType = NULL, chromosome.index = paste0("chr", c(1:22, "X")), species = "hg38") #species = "hg38", chromosome.index = paste0("chr", c(1:22, "X"))
# Add gene labels
circos.genomicLabels(plot_labels, labels.column = 4, side = "outside",  cex = 0.60, col = plot_labels$cols, line_col = plot_labels$cols)
# Add ideogram
circos.genomicIdeogram()
# Add Chr labels
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = gsub("chr", "", get.cell.meta.data("sector.index"))
  if (sector.name == "X") cols <- "cyan4" else if ((as.numeric(sector.name) %% 3) == 0) cols <- "cyan4" else if  ((as.numeric(sector.name) %% 2) == 0) cols <- "#0D0887FF" else cols <- "#A62098FF"
  circos.text(CELL_META$xcenter, 
              ylim[1] + cm_h(11), 
              sector.name, 
              facing = "clockwise",
              niceFacing = TRUE, 
              adj = c(0, 2),
              cex = 1,
              font = 2,
              col = cols)
}, bg.border = NA)
# Mutation Rate
circos.track(ylim = c(0, 90), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  
  y = seq(0, 90, by = 30)
  circos.segments(0, y, xlim[c("max.data")], y)
  
}, track.height=0.20, bg.border = NA)
circos.trackLines(track.index=4,final$chr, final$start, final$Rate, type="h", col = final$colour_rate)
# Ave Mutation Size
circos.track(ylim = c(-20, 40), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  
  y = seq(-20, 40, by = 10)
  circos.segments(0, y, xlim[c("max.data")], y, col="gray")
  
}, track.height=0.35, bg.border = NA)
circos.trackPoints(track.index=5, ave_size_2$chr, ave_size_2$start, ave_size_2$AVE, cex = 0.2, col = ave_size_2$col)#ave_size_2$Event_Type[ave_size_2$AVE > -20] )
dev.off()




