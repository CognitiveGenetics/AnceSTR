# Import Library
library(circlize)

# Import Data
# Import data
bed = generateRandomBed()

# Idiogram plot paraeters
circos.par("start.degree" = 90)
circos.initializeWithIdeogram(species = "hg38")

# Subset chromosomes
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(3,5,2,8)))

# Add track
circos.track(ylim = c(0, 100), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
}, track.height=0.3, bg.border = NA)

# Add plot
circos.trackLines(all_data$chr, all_data$start, all_data$Rate)
circos.trackPoints(bed$chr, bed$start, bed$value1*30)

# Reset plot
circos.clear()
