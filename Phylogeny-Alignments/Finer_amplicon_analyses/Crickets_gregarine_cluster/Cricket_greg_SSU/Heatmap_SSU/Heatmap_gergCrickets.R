# Load the necessary packages
library(gplots)
library(ggplot2)
library(extrafont)
library(viridis)
loadfonts()
setwd("C:/Users/eb826/OneDrive - University of Exeter/PhD_EdouardMicrosporidia/Thesis/Chapter3_Protist parasites survey/Crickets_greg")

# Read the CSV file into a data frame
df <- read.csv("blast_results.csv", row.names = 1)

# Print the data frame
print(df)

#windows(width = 10, height = 10)

png(filename="heatmap.png", res = 500, width = 10.2, height = 7, units = "in")

# For changing the font size of the legend title,labels and axis, the user needs to set
par(cex.main=1, cex.lab=1, cex.axis=1)

# Create a heatmap with dendrograms
heatmap.2(as.matrix(df), 
          trace = "none", 
          margins = c(17, 24),  # Reduce the size of the margins
          col = viridis(256),
          key=T,
          key.title = "Percentage Identity", 
          key.xlab = "Higher values indicate higher identity",
          #key.par =  list(cex=0.5),  # Reduce the margins around the key
          keysize=0.1,
          cexRow = 0.9, 
          cexCol = 0.9,
          offsetRow = 0.1,
          offsetCol = 0.1,
          density.info = "none", 
          lhei = c(1.7, 5.7), 
          lwid = c(1.4, 4.5),
          srtCol = 45)

#title(main="Percentage similiraties for gregarine amplicons from cricket hosts", line=0.3, adj=0.9) # Add a main title

# Stop the graphics device and save the image
dev.off()

png(filename="heatmap1.png", res = 500, width = 10.2, height = 7, units = "in")

# Create a heatmap without dendrograms
heatmap.2(as.matrix(df), 
          trace = "none", 
          margins = c(17, 24),  # Reduce the size of the margins
          col = viridis(256),
          key=T,
          key.title = "Percentage Identity", 
          key.xlab = "Higher values indicate higher identity",
          keysize=0.1,
          cexRow = 0.9, 
          cexCol = 0.9,
          offsetRow = 0.1,
          offsetCol = 0.1,
          density.info = "none", 
          lhei = c(1.7, 5), 
          lwid = c(1.4, 4.5),
          srtCol = 45,
          Rowv = FALSE,  # Disable the row dendrogram
          Colv = FALSE)  # Disable the column dendrogram

title(main="Percentage similiraties for gregarine amplicons from cricket hosts", line=0.5, adj=0.8) # Add a main title

# Stop the graphics device and save the image
dev.off()
