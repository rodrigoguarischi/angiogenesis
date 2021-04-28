### RNAseq-TLDA_comparison_plots.R ################################################################
# This script will make plots to compare results available from RNAseq with those obtained by
# Taqman Low Density Arrays (TLDA)

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd("~/angiogenesis_transcriptome/scripts/");

# Adds some localy installed packages to path
.libPaths( "/home/rodrigoguarischi/R/x86_64-pc-linux-gnu-library/3.3" );

# Load required libraries and datasets
library( BoutrosLab.plotting.general );
library( reshape2 );
library( org.Mm.eg.db );

### FUNCTIONS #####################################################################################

### Function - factor2numeric() #######################################################################
# Input variables: A factor with values to be converted
#
# Output variables: A vector with numeric values
#
# Description: Convert a factor to a vector of numeric values
factor2numeric <- function( x ) {
  
  stopifnot( inherits( x, "factor" ) );

  return( as.numeric( levels(x) )[x] );
  }

### DATA ANALYSIS #################################################################################

### Make scatter plot of genes by CTs

# Read TLDA results by well
tlda.wells <- read.csv(
  file = "aux_files/TLDA/Well Results.csv",
  header = TRUE, 
  comment.char = "#"
  );

# Read TLDA results file by biogroups
tlda.results <- read.csv(
  file = "aux_files/TLDA/Biogroup Results.csv",
  header = TRUE, 
  comment.char = "#"
  );

# Leave only Gene Symbol as target name 
tlda.wells$Target.Name <- gsub(
  "-(Mm|Hs).*$",
  "",
  tlda.wells$Target.Name
  );

tlda.results$Target.Name <- gsub(
  "-(Mm|Hs).*$",
  "",
  tlda.results$Target.Name
  );

# Filter wells that were ommited in the preprocessing steps
tlda.wells <- tlda.wells[tlda.wells$Omit == "false", ];

# Subset columns of interest
tlda.wells <- data.frame(
  Gene = tlda.wells$Target.Name,
  CT = factor2numeric( droplevels( tlda.wells$Eq..Cq ) ),
  Condition = tlda.wells$Biological.Group.Name,
  Sample = factor(tlda.wells$Sample.Name)
  );

# Subset TLDA results to columns of interest and cast to appropriated data type (Warnings to Endogenous Control genes will appear, ok!)
tlda.results <- data.frame(
  Condition = factor( tlda.results$Bio.Group.Name ),
  Gene = factor( tlda.results$Target.Name ),
  Rq = factor2numeric( droplevels(tlda.results$Rq) ),
  Mean.CT = tlda.results$Mean.Equivalent.Cq
  );

# Transform Rq values to log2 scale
tlda.results$Rq <- log2(tlda.results$Rq);

# # Transform data into an easier data.frame (Genes by Samples)
# acast(
#   data = tlda.results, 
#   formula = Gene ~ Condition,
#   value.var = "Rq"
#   );

# Create object to gather all graph parameters
graph.parameters <- list();

# Sort genes by CT values
graph.parameters$CT.ordered <- names(
  sort(
    tapply(
      X = tlda.wells$CT, 
      INDEX = tlda.wells$Gene,
      FUN = max
      )
    )
  );

# Determine the xaxis position for each gene sorting by maximum CT value
graph.parameters$xpos <- 1:length( unique( tlda.wells$Gene ) );
names(graph.parameters$xpos) <- graph.parameters$CT.ordered;

# Mark genes that have P12 (reference sample) with mean CT above threshold as group 2
graph.parameters$cycle.threshold <- 35;
graph.parameters$groups <- rep( 1, length( graph.parameters$CT.ordered ) );
graph.parameters$groups[ graph.parameters$CT.ordered %in% tlda.results[ tlda.results$Condition == 'P12' & tlda.results$Mean.CT > graph.parameters$cycle.threshold,]$Gene ] <- 2
names(graph.parameters$groups) <- graph.parameters$CT.ordered;

# Add column Xpos and Groups to tlda.wells
tlda.wells$Xpos <- graph.parameters$xpos[ as.character( tlda.wells$Gene ) ];

# Make scatterplot of genes by CT values sorted by CT value for the gene
CTs.scatterplot <- create.scatterplot(
  #
  formula = CT ~ Xpos,
  data = tlda.wells,
  # main = 'Genes by CTs in qRT-PCR',
  # main.cex = 2,

  # Define point style and add colours to plot
  pch = 1,
  col = 'black',

  # Configure x-axis
  xlab.label = "Genes",
  # xlab.cex = 3.5,
  xlimits = c(1, max(tlda.wells$Xpos)) + c(-1,1),
  xat = 1:max(tlda.wells$Xpos),
  xaxis.cex = 1.3,
  xaxis.fontface = 1,
  xaxis.lab = graph.parameters$CT.ordered,
  xaxis.rot = 90,
  
  # Configure y-axis
  ylab.label = "CT",
  # ylab.cex = 3.5,
  yaxis.cex = 1.3,
  yaxis.fontface = 1,

  # # Add horizontal line
  # abline.h = graph.parameters$cycle.threshold,
  # abline.col = "red",
  # abline.lwd = 3,
  # abline.lty = 2,

  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# tiff(
#   filename = "aux_files/plots/tlda_cts.tiff",
#   type = "cairo",
#   width = 30,
#   height = 15,
#   units = 'cm',
#   res = 500,
#   compression = 'lzw',
#   bg = NA
#   );
# CTs.scatterplot;
# dev.off();

# Make plot of observed diferences (in CTs) between pipeting replicates
experiment.replicates <- factor( paste(tlda.wells$Gene, tlda.wells$Sample, sep = "-") );
replicates.deltaCTs <- sapply( 
  split( 
    seq( along = experiment.replicates ),
    experiment.replicates
    ),
  function(rows) if( length(rows) > 1 ) diff( tlda.wells[rows,"CT"] ) else 0
  );

deltaCTs.plot.data <- data.frame(
  deltaCTs = abs(replicates.deltaCTs),
  Gene = gsub(
    "-.*$",
    "",
    names( replicates.deltaCTs )
    ),
  row.names = NULL
  );

deltaCTs.plot.data$Xpos <- graph.parameters$xpos[ as.character( deltaCTs.plot.data$Gene ) ];

replicates.deltaCTs.scatterplot <- create.scatterplot(
  #
  formula = deltaCTs ~ Xpos,
  data = deltaCTs.plot.data,

  # Define point style and add colours to plot
  pch = 21,
  col = 'black',
  
  # Configure x-axis
  xlab.label = "Genes",
  xlab.cex = 1.5,
  xlimits = c(1, max(tlda.wells$Xpos)) + c(-1,1),
  xat = 1:max(tlda.wells$Xpos),
  xaxis.cex = 0.7,
  xaxis.fontface = 1,
  xaxis.lab = rep(' ', max(tlda.wells$Xpos)),

  # Configure y-axis
  ylab.label = expression( paste( "Replicates ", Delta, "CT") ),
  ylab.cex = 1.5,
  yaxis.cex = 1,
  yaxis.fontface = 1,
  ylimits = c(0,4),
  yat = seq(0,4,1),

  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Aggregate qRT-PCR plots into one multiplot
CTs.multiplots <- create.multiplot(
  #
  plot.objects = list(CTs.scatterplot, replicates.deltaCTs.scatterplot),
  panel.heights = c(0.25,0.75),

  # Configure x-axis
  xlab.label = c("Genes"),
  xlab.cex = 1,
  xaxis.rot = 90,
  xaxis.fontface = 1,
  xaxis.cex = 0.7,
  
  # Configure y-axis
  ylab.label = c( expression( paste( "Replicates ", Delta, "CT") ), "\t" , "CT", "\t"),
  ylab.cex = 1.1,
  yaxis.fontface = 1,
  yaxis.cex = 0.7,
  y.spacing = -2,
  retrieve.plot.labels = TRUE,
  
  description = 'Multiplot created by BoutrosLab.plotting.general'
  );

### Export scatterplots into a png file
png(
  filename="~/Desktop/tlda-results.CT35.png",
  type="cairo",
  units="in",
  width = 10,
  height = 8,
  pointsize=12,
  res=500
  );

plot( CTs.multiplots );

dev.off();


### Make scatter plot of fold-change in RNAseq and TLDA to compare them

# Call build_ddsFull.R script to get a DESeqDataSet object
source( file = "build_ddsFull.R" );

# Create factor that will be used in contrast
dds <- ddsFull;
dds$Group <- factor( paste0(dds$Treat, dds$Time) );
design(dds) <- ~ Group;
dds <- DESeq(dds);

# Remove list of genes not evaluated by RNAseq and endogenous control genes (rRNA 18S, Sdha and Tbp).
# Muc2 gene was splitted in two on new ENSEMBL annotation so we decided to remove it from analysis as well!
graph.parameters$fold.change.genes <- graph.parameters$CT.ordered[ -which( graph.parameters$CT.ordered %in% c("18S","Rn18s", "Sdha", "Tbp", "Muc2") ) ]; 

# Build data.frame for RNAseq results
rnaseq.results <- sapply( 
  c("P15", "P17", "R12", "R12.5", "R15", "R17"),
  function(sample.name) results(
    dds, 
    contrast = c("Group", sample.name, "P12")
    )[ mapIds( x = org.Mm.eg.db, keys = graph.parameters$fold.change.genes, keytype = "SYMBOL", column = "ENSEMBL", multiVals="first" ), "log2FoldChange"]
  );
rownames(rnaseq.results) <- graph.parameters$fold.change.genes;

# Melt the data frame into angio.results and add TLDA results afterwards. 
# Work with only genes in gene.list and remove reference condition (P12, Rq=1)
angio.results <- melt( rnaseq.results );
colnames(angio.results) <- c("Gene", "Condition", "RNASeq");
angio.results$TLDA <- apply(
  angio.results, 
  1, 
  function(df.row) subset( tlda.results, tlda.results$Gene == df.row["Gene"] & tlda.results$Condition == df.row["Condition"] )$Rq
  );

# Calculate scatter plot boundaries
graph.parameters$axis.limits <- floor( range( angio.results[, c("TLDA","RNASeq") ] ) ) + c(0,1);

# Get list of experiments that are above threshold
graph.parameters$high.CT.experiments <- as.character(
  apply(
    tlda.wells[ tlda.wells$CT > graph.parameters$cycle.threshold , c("Gene", "Condition")],
    1,
    paste,
    collapse = '-'
    )
  );
  
# Assign a group to each experiment. 1 = bellow CT threshold and 2 = above CT threshold
angio.results$Group <- graph.parameters$groups[ as.character(angio.results$Gene) ];
angio.results[ paste( angio.results$Gene, angio.results$Condition, sep = "-") %in% graph.parameters$high.CT.experiments, ]$Group <- 2;

# Generate scatterplot for experiments that are bellow CT threshold
fold.change.scatterplot.g1 <- create.scatterplot(
  #
  formula = TLDA ~ RNASeq,
  data = subset(angio.results, angio.results$Group == 1),
  type = c("p", "r"),

  # Customize points style  
  pch = 21,
  col = "black",
  
  # add x=y line
  add.xyline = TRUE,
  xyline.lwd = 2,
  xyline.lty = 2,
  xyline.col = 'seagreen4',
  
  # Draw x = 0 and y = 0 lines
  abline.h = 0,
  abline.v = 0,
  abline.lwd = 1,
  
  # Change x axis
  xlimits = graph.parameters$axis.limits,
  xlab.label = "RNASeq (log2 fold-change)",
  xlab.cex = 1.3,
  xaxis.cex = 1,
  xaxis.fontface = 1,
  xaxis.rot = 0,
  xaxis.tck = c(1,-1),
  
  # Change y axis
  ylimits = graph.parameters$axis.limits,
  ylab.label = "qRT-PCR (log2 fold-change)",
  ylab.cex = 1.3,
  yaxis.cex = 1,
  yaxis.fontface = 1,
  yaxis.tck = c(1,-1),

  # Adding legend
  print.new.legend = TRUE,
  legend = list (

    # Adding correlation key to plot 1
    inside = list(
      fun = draw.key,
      args = list(
        key = get.corr.key(
          x = subset(angio.results, angio.results$Group == 1)$RNASeq,
          y = subset(angio.results, angio.results$Group == 1)$TLDA,
          label.items = c('spearman', "spearman.p"),
          alpha.background = 0
          )
        ),
      x = 0.02,
      y = 0.97,
      corner = c(0,1)
      )
    ),
  
  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

### Export scatterplot of "good" experiments into a png file

png(
  filename="~/Desktop/tlda-rnaseq_correlation.CT35.png",
  type="cairo",
  units="in",
  width = 8,
  height = 8,
  pointsize=12,
  res=500
  );

plot( fold.change.scatterplot.g1 );

dev.off();

# Generate scatterplot for experiments that are above CT threshold
fold.change.scatterplot.g2 <- create.scatterplot(
  #
  formula = TLDA ~ RNASeq,
  data = subset(angio.results, angio.results$Group == 2),
  type = c("p", "r"),
  main = 'qRT-PCR/RNAseq comparison',
  main.cex = 2,
  
  # Customize points style  
  pch = 21,
  col = graph.parameters$colour.pallet[2],
  
  # add x=y line
  add.xyline = TRUE,
  xyline.lwd = 2,
  xyline.lty = 2,
  xyline.col = 'seagreen4',
  
  # Draw x = 0 and y = 0 lines
  abline.h = 0,
  abline.v = 0,
  abline.lwd = 1,
  
  # Change x axis
  xlimits = graph.parameters$axis.limits,
  xlab.label = "RNASeq (log2 fold-change)",
  xlab.cex = 1.3,
  xaxis.cex = 1,
  xaxis.fontface = 1,
  xaxis.rot = 0,
  xaxis.tck = c(1,-1),
  
  # Change y axis
  ylimits = graph.parameters$axis.limits,
  ylab.label = "qRT-PCR (log2 fold-change)",
  ylab.cex = 1.3,
  yaxis.cex = 1,
  yaxis.fontface = 1,
  yaxis.tck = c(1,-1),
  
  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Create legend to be used in plot
legend.grob <- legend.grob(
  legends = list(
    legend = list(
      colours = graph.parameters$colour.pallet,
      labels = paste0( c('< ', '> '), graph.parameters$cycle.threshold, ' cycles' ),
      title = 'CT threshold',
      title.just = 'left',
      border = 'white'
    )
  ), 
  size = 1.5
  );

# Create multiplot
fold.change.multiplot <- create.multiplot(
  #
  plot.objects = list(fold.change.scatterplot.g2, fold.change.scatterplot.g1),
  panel.heights = c(0.5,0.5),
  main = "qRT-PCR/RNAseq comparison",
  main.cex = 2,
  
  # Configure x-axis
  xlab.label = c("RNASeq (log2 fold-change)"),
  xlab.cex = 1,
  xaxis.fontface = 1,
  xaxis.cex = 0.7,
  
  # Configure y-axis
  ylab.label = c( "qRT-PCR (log2 fold-change)", "qRT-PCR (log2 fold-change)"),
  ylab.cex = 1,
  yaxis.fontface = 1,
  yaxis.cex = 0.7,
  retrieve.plot.labels = TRUE,
  
  # Adding legend
  print.new.legend = TRUE,
  legend = list (
    right = list( fun = legend.grob ),
    
    # Adding correlation key to plot 1
    inside = list(
      fun = draw.key,
      args = list(
        key = get.corr.key(
          x = subset(angio.results, angio.results$Group == 1)$RNASeq,
          y = subset(angio.results, angio.results$Group == 1)$TLDA,
          label.items = c('pearson','pearson.p', 'spearman', "spearman.p"),
          alpha.background = 0,
          key.cex = 0.75
          )
        ),
        x = 0.02,
        y = 0.97,
        corner = c(0,1)
      ),
    
    # Adding correlation key to plot 2
    inside = list(
      fun = draw.key,
      args = list(
        key = get.corr.key(
          x = subset(angio.results, angio.results$Group == 2)$RNASeq,
          y = subset(angio.results, angio.results$Group == 2)$TLDA,
          label.items = c('pearson','pearson.p', 'spearman', "spearman.p"),
          alpha.background = 0,
          key.cex = 0.75
          )
        ),
        x = 0.02,
        y = 0.43,
        corner = c(0,1)
      )
    ),

  description = 'Multiplot created by BoutrosLab.plotting.general'
  );


### Export scatterplots into a png file side by side
png(
  filename="~/Desktop/tlda-rmaseq_comparison.CT35.png",
  type="cairo",
  units="in",
  width = 16,
  height = 8,
  pointsize=12,
  res=500
  );

plot(
  CTs.multiplots,
  position = c(0,0,0.56,0.995),
  newpage = TRUE
  );

print(
  fold.change.multiplot,
  position = c(0.57,0.044,1,1), 
  newpage = FALSE
  );

dev.off();

