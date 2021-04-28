### literature_signatures_comparison.R ############################################################
# This script is used to make comparisons between existing literature signatures (angiogenic and 
# hypoxic). It takes some tabular files as input and generate a series of plots.

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd("~/angiogenesis_transcriptome/scripts/");

# Load required libraries and datasets
library( VennDiagram );
library( org.Hs.eg.db );
library( gtools );
library( reshape2 );
library( gridExtra );
library( limma );       # Only if standarizeGeneIDs() is used
library( biomaRt );     # Only if standarizeGeneIDs() is used
library( scales );

# Loading BL internal libraries (R-BL as of April 5, 2017)
library( BoutrosLab.statistics.survival );
library( BoutrosLab.datasets.breast.cancer );
library( BoutrosLab.plotting.general );
library( BoutrosLab.plotting.survival );

### Function - returnOldIfInvalid() ###############################################################
# Description: Function that check if new value is valid, otherwise returns old value.
returnOldIfInvalid <- function (value.old, value.new){
  
  if( length(value.new) == 1 ){
    if( ! is.na(value.new) ){
      return( as.character(value.new) );  
      }
    else{
      return( as.character(value.old) );
      }
    }
  else{
    return( as.character(value.old) );
   }
  }

### Function - standarizeGeneIDs() ################################################################
# Input variables: data.frame   First column should have geneIDs to be verified
#
# Output variables: data.frame   Rownames will receive oficial geneIDs

# Description: Function that gets oficial geneIDs from a list of geneIDs. If geneIDs are from 
# Mus musculus get its respective human homologous from biomaRt package
standarizeGeneIDs <- function( input.df, organism="Hs" ){
  
  stopifnot( inherits( input.df, "data.frame" ) );
  
  switch(
    organism,
    Hs={
      rownames( input.df ) <- lapply( 
        input.df[,1], 
        function(x) returnOldIfInvalid( 
          x, 
          alias2Symbol(x, species = "Hs")
          )
        );
      },
    Mm={
      input.df$OficialMouseSymbol <- lapply( 
        input.df[,1], 
        function(x) returnOldIfInvalid(
          x,
          alias2Symbol(x, species = "Mm")
          )
        );
      
      human.mart <- useMart( biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org" );
      mouse.mart <- useMart( biomart="ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" );
      
      homologous <- getLDS(
        attributes = "mgi_symbol",
        filters = "mgi_symbol",
        values = input.df$OficialMouseSymbol,
        mart = mouse.mart,
        attributesL = "hgnc_symbol",
        martL = human.mart
        );
      
      mmusculus.genes <- factor( homologous$MGI.symbol );
      mmusculus2hsapiens <- sapply( 
        split( 
          seq( along = mmusculus.genes ),
          mmusculus.genes
          ),
        function(rows) paste( homologous[rows, ]$HGNC.symbol, collapse = '-' )
        );
      
      rownames( input.df ) <- sapply(
        input.df$OficialMouseSymbol, 
        function(x) returnOldIfInvalid(
          x,
          mmusculus2hsapiens[x]
          )
        );

      },
    {
      stop("Invalid organism. Possible values are \"Hs\" or \"Mm\".")      
      }
    )
  
  return( input.df );
  }

### Function - convertIDs() #######################################################################
# Input variables: List of IDs from one source (e.g: entrezIDs, symbol, Ensembl)
#
# Output variables: List of converted IDs
#
# Description: Convert a list of IDs from one identifier to another using AnnotationDbi package
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  
  stopifnot( inherits( db, "AnnotationDb" ) );
  ifMultiple <- match.arg( ifMultiple );
  
  suppressWarnings( selRes <- AnnotationDbi::select( db, keys=ids, keytype=fromKey, cols=c(fromKey,toKey) ) );
  
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ];
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ];
    }
  
  return( selRes[ match( ids, selRes[,1] ), 2 ] );
  }

### Function - mergeDataframesByRownames() ########################################################
# Input variables: Two data frames
#
# Output variables: One merged data frame, by rownames
#
# Description: Merge two dataframes by rowname, replace its rownames and remove duplicated column
mergeDataframesByRownames <- function( dataframe.A, dataframe.B ){

  stopifnot( inherits( dataframe.A, "data.frame" ) );
  stopifnot( inherits( dataframe.B, "data.frame" ) );
  
  merged.dataframe <- merge(
    dataframe.A,
    dataframe.B,
    by = "row.names",
    all = TRUE
    );
  
  rownames(merged.dataframe) <- merged.dataframe$Row.names;
  merged.dataframe <- merged.dataframe[, -which(colnames(merged.dataframe) == "Row.names"), drop=FALSE];
  
  return(merged.dataframe);
  }

### Function - calcHazardRatios() #################################################################
# Input variables: Survival object, respective sample group for each sample and signature's name
#
# Output variables: One data frame with HR, lower_95CI, upper_95CI, pval, nevents
#
# Description: Calculate Hazard Ratios from a survival object and respective sample group
calcHazardRatios <- function( surv.obj, sample.group, signature.name ){
  
  stopifnot( inherits( surv.obj, "Surv" ) );
  stopifnot( ! length( unique( sample.group ) ) != 2 ); # Confirms that only 2 classes are given
  stopifnot( inherits( signature.name, "character" ) );
  
  coxph <- coxph( surv.obj ~ sample.group );
  
  coxph.summary <- summary( coxph );
  
  if( ph.fails( cox.model = coxph, cox.zph.threshold = 0.1, pvalues = FALSE ) == TRUE ) {
    warning( paste0("PH hazard assumption fails for signature: ", signature.name) );
  };
  
  hazard.ratios <- data.frame(
    HR = coxph.summary$conf.int[,1],
    lower_95CI = coxph.summary$conf.int[,3],
    upper_95CI = coxph.summary$conf.int[,4],
    # pval = 2 * pnorm( -abs( coxph.summary$coefficients[,4] ) ),
    pval = coxph.summary$coefficients[,5],
    nevent = coxph.summary$nevent,
    row.names = signature.name
    );
  
  return(hazard.ratios);
  }

### DATA ANALYSIS #################################################################################

### Prepare a table of genes by signatures

# Initialize object 'signatures' that will hold data from all signatures
signatures <- list();

# Initialize empty data.frame to be populated
signatures$data <- data.frame();

# Initialize empty array. This array will be used to discriminate angiogenic (A) and hypoxic (H) signatures columns
signatures$type <- character();

# Initialize object 'tmp' to hold all temporary objects
tmp <- list();

# Get list of signatures' files
tmp$signatures.files <- c(
  list.files( 
    path="aux_files/angiogenesis_signatures", 
    full.names = TRUE, 
    pattern="*\\.txt$" ), 
  list.files( 
    path="aux_files/hypoxia_signatures", 
    full.names = TRUE, 
    pattern="*\\.txt$" )
  );

# Load signatures into data frame
for( file in tmp$signatures.files ){
  
  # Get signature's name (from filename)
  tmp$signature.name <- gsub( "\\.txt", "", basename(file) );
  
  # Transform special character encoded in filenames (SVN has limitations with special characters)
  tmp$signature.name <- gsub( "-open.bracket-", "\\[", tmp$signature.name );
  tmp$signature.name <- gsub( "-close.bracket-", "\\]", tmp$signature.name );
  tmp$signature.name <- gsub( "-percent-", "\\%", tmp$signature.name );

  # Get signature's type (from file path)
  if( grepl( "angiogenesis_signatures", file ) ) { signatures$type[[ tmp$signature.name ]] <- "A" }
  else if( grepl( "hypoxia_signatures", file ) ) { signatures$type[[ tmp$signature.name ]] <- "H" } 
  else { warning( paste0("Invalid signature type for file: ", file) ) };
  
  # Load into temporary object
  tmp$signature <- read.table (
    file, 
    sep = "\t", 
    header = TRUE, 
    comment.char = "#"
    );
  
  # Transform data.frame the way standarizeGeneIDs function expects it
  if( "A" == signatures$type[[tmp$signature.name]] ) {
    tmp$signature <- tmp$signature[ , c("InputMatch", "Expression") ];
  } else if( "H" == signatures$type[[tmp$signature.name]] ) {
    tmp$signature <- data.frame (
      InputMatch = tmp$signature$Symbol,
      Expression = TRUE
      );
  } 
  
  # Segregate signature's origin (Human or Mouse)
  if( tmp$signature.name != "Langlois_2014" ){
    tmp$signature <- standarizeGeneIDs( tmp$signature, organism = 'Hs' )[,"Expression", drop = FALSE];
    } else {
    tmp$signature <- standarizeGeneIDs( tmp$signature, organism = 'Mm' )[,"Expression", drop = FALSE];
    }
  
  # Assign signature's name to column
  colnames(tmp$signature) <- tmp$signature.name;
  
  # Combine 
  signatures$data <- mergeDataframesByRownames( signatures$data, tmp$signature );
  }

# Get number of genes in each signature
signatures$n.genes <- colSums( ! is.na(signatures$data) );

# jaccard.distance( t( ! is.na(signatures$data) ) );
# dim(signatures$data);

### Create a barplot of gene symbols by number of signatures

# Initialize object 'barplot' that will hold data for barplot
barplot <- list();

# Build data.frame to use in barplot with data from signatures$data (angiogenesis-only)
barplot$angio.sigs.gene.frequency <- rowSums( ! is.na( signatures$data[, signatures$type == 'A'] ) );
barplot$angio.sigs.gene.frequency <- barplot$angio.sigs.gene.frequency[ barplot$angio.sigs.gene.frequency >= 1 ];
signatures$angiogenesis.geneSymbols <- names( barplot$angio.sigs.gene.frequency );

# Subset genes within at least 2 different signatures 
barplot$data <- data.frame(
  Symbol = names( barplot$angio.sigs.gene.frequency[ barplot$angio.sigs.gene.frequency >= 2 ] ),
  Frequency = barplot$angio.sigs.gene.frequency[ barplot$angio.sigs.gene.frequency >= 2 ],
  stringsAsFactors = FALSE
  );

# Sort barplot using Frequency as decreasing and Symbol as ascending keys
barplot$data <- barplot$data[ order( barplot$data$Frequency, rev(barplot$data$Symbol), decreasing = TRUE), ];
rownames(barplot$data) <- NULL;

# Get list of genes in hypoxia signature to colour barplot
barplot$hypoxia.sigs.gene.frequency <- rowSums( ! is.na( signatures$data[, signatures$type == 'H'] ) );
barplot$hypoxia.sigs.gene.frequency <- barplot$hypoxia.sigs.gene.frequency[ barplot$hypoxia.sigs.gene.frequency >= 1 ];
signatures$hypoxia.geneSymbols <- names( barplot$hypoxia.sigs.gene.frequency );

# Set colours for bars
barplot$pallet <- c('black', 'darkgray');
barplot$bar.colours <- rep( barplot$pallet[1], nrow(barplot$data) );
barplot$bar.colours[ barplot$data$Symbol %in% signatures$hypoxia.geneSymbols ] <- barplot$pallet[2];

# Uses a trick to rotate plot without messing up data
barplot$plot <- create.barplot(
  #
  formula =  rev( as.numeric( rownames(barplot$data) ) ) ~ Frequency,
  data = barplot$data,

  plot.horizontal = TRUE,
  col = rev( barplot$bar.colours ),
  border.col = 'transparent',
  
  left.padding = 5,
  
  xlab.lab = 'Number of signatures',
  xaxis.fontface = 1,
  xlimits = c(0, max(barplot$data$Frequency)*1.1),
  xat = seq(0, max(barplot$data$Frequency)*1.1, 1),

  ylab.lab = 'Gene Symbol',
  yaxis.fontface = 1,
  yaxis.lab = rev( barplot$data$Symbol ),

  # Adding legend to explain bar colour-coding
  legend = list(
    right = list( 
      fun = legend.grob(
        legends = list(
          legend = list(
            colours = barplot$pallet,
            labels = c('Not present', 'Also present'),
            title = " Hypoxia signatures",
            border = 'white'
            )
          ),
          title.just = 'left',
          title.cex = 1.25,
          label.cex = 1.15
        )
      )
    ),

  description = 'Barplot created by BoutrosLab.plotting.general'
  );

# Export image to tiff
tiff(
  filename = "aux_files/plots/signatures_meta_analysis_a.tiff", 
  type = "cairo", 
  width = 20, 
  height = 25, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

plot(barplot$plot);

dev.off();


### Create a venn diagram to compare angiogenesis signatures with hypoxia ones

venn.plot <- venn.diagram(

  x = list(
    hypoxia = signatures$hypoxia.geneSymbols,
    angiogenesis = signatures$angiogenesis.geneSymbols
    ),
  filename = NULL,

  ext.percent = 0,
  lty = "blank",
  fill = c("dodgerblue","red"),
  alpha = 0.35,
  label.col = c("darkblue", "white", "darkred"),
  cex = 2,
  fontface = "bold",

  cat.col = c("darkblue", "darkred"),
  cat.cex = 2,
  cat.pos = c(15,-15),
  cat.default.pos = "text",
  cat.dist = 0.05
  );

# Remove log file created
file.remove( list.files( pattern="Venn.*.log$") )

#
# Perform Hypergeometric test to access gene overlap significancy
# Assuming number of genes equals to 20296 (# of coding genes in GRCh38.p3 at ENSEMBL statistics)
#

phyper(
  q = sum(signatures$angiogenesis.geneSymbols %in% signatures$hypoxia.geneSymbols), # Overlap between 2 sets
  m = length(signatures$angiogenesis.geneSymbols),                                  # Number of angiogenesis genes
  n = 20296 - length(signatures$angiogenesis.geneSymbols),                          # Total number of genes - number of angiogenesis genes
  k = length(signatures$hypoxia.geneSymbols),                                       # Number of hypoxia genes
  lower.tail = FALSE                                                                # What I want is P[X > x].
  );

#
# p = 8.945098e-05
#

# Export image to tiff
tiff(
  filename = "aux_files/plots/signatures_meta_analysis_c.tiff", 
  type = "cairo", 
  width = 15, 
  height = 15, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

# # Display plot
# grid.newpage();
grid.draw(venn.plot);

dev.off();

### Create a heatmap to compare angiogenesis signatures with hypoxia ones

# Initialize object 'shared.genes' that will hold data for signatures' heatmap of shared genes
shared.genes <- list(); 

# Calculate all possible permutation of 2 columns
shared.genes$data <- as.data.frame( 
  permutations( 
    ncol(signatures$data),
    2,
    repeats.allowed = TRUE
    )
  );

# Replace column index by signature's names
shared.genes$data$V1 <- colnames(signatures$data)[shared.genes$data$V1];
shared.genes$data$V2 <- colnames(signatures$data)[shared.genes$data$V2];

# Add column to data frame with number of shared genes between each pair of columns
shared.genes$data$shared.genes <- apply( 
  shared.genes$data, 
  1, 
  function(signatures.pairs) sum( rowSums( ! is.na( signatures$data[, signatures.pairs] ) ) == 2 )
  );

# Reshape dataframe object
shared.genes$data <- acast(
  data = shared.genes$data,
  formula = V2 ~ V1,
  value.var = "shared.genes"
  );

# Calculate the respective percentage it correspond from the number of genes in signatures
shared.genes$data <- shared.genes$data / ( signatures$n.genes[ colnames( shared.genes$data ) ] %*% t( rep(1, ncol(shared.genes$data) ) ) );

# Define a colours to be used in heatmap (use same colour scheme from venn diagram on covariates)
shared.genes$cells.pallet <- c('white','firebrick');
shared.genes$covariates.pallet <- adjustcolor( c("red", "dodgerblue"), 0.35 );

# Make covariates
shared.genes$cov.colours <- rep( shared.genes$covariates.pallet[1], nrow(shared.genes$data) );
shared.genes$cov.colours[ colnames( shared.genes$data ) %in% names( signatures$type["H"  == signatures$type ] ) ] <- shared.genes$covariates.pallet[2];

shared.genes$covariates <- list(
  rect = list(
    col = 'white',
    fill = shared.genes$cov.colours,
    lwd = 1.5
    )
  );

# join covariate legends
shared.genes$legend <- list(
  legend = list(
    colours = shared.genes$covariates.pallet,
    labels = c('angiogenesis','hypoxia'),
    title = 'Signature type',
    border = 'white'
  ),
  legend = list(
    colours = shared.genes$cells.pallet,
    labels = c('0%','100%'),
    title = 'Shared genes (%)',
    continuous = TRUE,
    height = 5,
    pos.x = 0.42
    )
  );

# Dendrogram provided
shared.genes$dendrogram <- create.dendrogram(
  x = shared.genes$data,
  cluster.dimension = 'row'
  );

# Plot shared genes heatmap
shared.genes$heatmap <- create.heatmap(
  #
  x = shared.genes$data,
  colour.scheme = shared.genes$cells.pallet,
  
  # Use same clustering for rows and columns
  clustering.method = 'none',
  row.dendrogram = shared.genes$dendrogram,
  col.dendrogram = shared.genes$dendrogram,

  top.dendrogram.size = 1.25,
  right.dendrogram.size = 1.25,
  
  # Active right-side and top covariate
  covariates = shared.genes$covariates,
  covariates.top = shared.genes$covariates,

  # Disbale colour.key and print legend
  print.colour.key = FALSE,
  covariate.legend = shared.genes$legend,
  legend.side = 'right',
  legend.cex = 0.5,
  legend.title.just = 'left',
  legend.title.cex = 0.5,

  axes.lwd = 1,
  
  xaxis.lab = gsub("_", " ", colnames( shared.genes$data ) ),
  xaxis.fontface = 1,
  xaxis.cex = 0.5,
  
  yaxis.lab = gsub("_", " ", colnames( shared.genes$data ) ),
  yaxis.fontface = 1,
  yaxis.cex = 0.5,

  description = 'Heatmap created using BoutrosLab.plotting.general'
  );

# Export image to tiff
tiff(
  filename = "aux_files/plots/signatures_meta_analysis_b.tiff", 
  type = "cairo", 
  width = 13, 
  height = 10, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

plot(shared.genes$heatmap);

dev.off();

### Create a heatmap to compare angiogenesis signatures with hypoxia ones

# Initialize object 'signatures.correlation' that will hold data for signatures' correlation heatmap
signatures.correlation <- list(); 

# Calculate correlation between all signatures (DISREGARD EXPRESSION INFORMATION!)
signatures.correlation$data <- cor( ! is.na(signatures$data) );

# Define a colours to be used in heatmap
signatures.correlation$cells.pallet <- c( 'red', 'white', 'turquoise' );
signatures.correlation$covariates.pallet <- adjustcolor( c('red', 'dodgerblue'), 0.35);

# Make covariates
signatures.correlation$cov.colours <- rep( signatures.correlation$covariates.pallet[1], nrow(signatures.correlation$data) );
signatures.correlation$cov.colours[ colnames(signatures.correlation$data) %in% names( signatures$type["H"  == signatures$type] ) ] <- signatures.correlation$covariates.pallet[2];

signatures.correlation$covariates <- list(
  rect = list(
    col = 'white',
    fill = signatures.correlation$cov.colours,
    lwd = 1.5
    )
  );

# join covariate legends
signatures.correlation$legend <- list(
  legend = list(
    colours = signatures.correlation$covariates.pallet,
    labels = c('angiogenesis','hypoxia'),
    title = 'Signature type',
    border = 'white'
    ),
  legend = list(
    colours = signatures.correlation$cells.pallet,
    labels = c('-1', '1'),
    title = 'Correlation',
    continuous = TRUE,
    height = 5,
    pos.x = 0.3
    )
  );

# Dendrogram provided
signatures.correlation$dendrogram <- create.dendrogram(
  x = signatures.correlation$data,
  cluster.dimension = 'row'
  );

# Plot correlation heatmap
signatures.correlation$heatmap <- create.heatmap(
  #
  x = signatures.correlation$data,
  colour.scheme = signatures.correlation$cells.pallet,
  colour.centering.value = 0,
  at = seq(-1, 1, 0.001),
  
  # Use same clustering for rows and columns
  clustering.method = 'none',
  row.dendrogram = signatures.correlation$dendrogram,
  col.dendrogram = signatures.correlation$dendrogram,

  top.dendrogram.size = 1.25,
  right.dendrogram.size = 1.25,
  
  # Active right-side and top covariate
  covariates = signatures.correlation$covariates,
  covariates.top = signatures.correlation$covariates,
  
  # Disbale colour.key and print legend
  print.colour.key = FALSE,
  covariate.legend = signatures.correlation$legend,
  legend.side = 'right',
  legend.cex = 0.5,
  legend.title.just = 'left',
  legend.title.cex = 0.5,
  
  axes.lwd = 1,
  
  xaxis.lab = gsub("_", " ", colnames(signatures.correlation$data) ),
  xaxis.fontface = 1,
  xaxis.cex = 0.5,
  
  yaxis.lab = gsub("_", " ", colnames(signatures.correlation$data) ),
  yaxis.fontface = 1,
  yaxis.cex = 0.5,
  
  description = 'Heatmap created using BoutrosLab.plotting.general'
  );

# Export image to tiff
tiff(
  filename = "aux_files/plots/signatures_meta_analysis_d.tiff", 
  type = "cairo", 
  width = 13, 
  height = 10, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

plot(signatures.correlation$heatmap);

dev.off();

### Create a multiplot to compare angiogenesis signatures using Metabric dataset
metabric <- load.breast.cancer.datasets( datasets.to.load=c('Metabric.Training','Metabric.Validation') );
# metabric$all.data$Metabric.Training[1:10,1:10]

# Initialize object 'metabric.analysis' that will hold data for signatures evaluation on Metabric dataset
metabric.analysis <- list(); 

# Combine two survival objects from metabric dataset (training/validation) into one
metabric.analysis$survival <- rbind( metabric$all.survobj$Metabric.Training, metabric$all.survobj$Metabric.Validation );
metabric.analysis$surv.obj <- Surv( metabric.analysis$survival[,1], metabric.analysis$survival[,2] );

# Truncate survival object at 15 years
metabric.analysis$surv.obj <- survival.truncate( list( metabric.analysis$surv.obj ), 15)[[1]];
metabric.analysis$survival <- metabric.analysis$surv.obj[ , c(1,2) ];

# Get percentile of events in dataset
metabric.analysis$percentile.of.events <- prop.table( table( metabric.analysis$surv.obj[,2] ) )[1];

# Retrieve mapping of symbols to geneIDs from BioMart
tmp$symbol2entrez <- getBM(
  attributes = c("hgnc_symbol", "entrezgene"),
  filters = "hgnc_symbol",
  values = rownames(signatures$data),
  mart = useMart( biomart="ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl" )
  );

# Initialize object to store sample scores for each signature
metabric.analysis$sample.scores <- data.frame(
  row.names = c( colnames( metabric$all.data$Metabric.Training ), colnames( metabric$all.data$Metabric.Validation ) )
  );

# Initialize object to store survival information for each signature
metabric.analysis$sign.odds.ratios <- data.frame(
  HR = numeric,
  lower_95CI = numeric,
  upper_95CI = numeric,
  pval = numeric,
  nevent = numeric
  );

# For each signature calculate sample's score and make survival analysis
for ( signature.name in colnames( signatures$data[, signatures$type == 'A' ] ) ){

  # Get list of genes present in the signature
  tmp$signature <- subset(tmp$symbol2entrez, hgnc_symbol %in% rownames( signatures$data[ ! is.na( signatures$data[[signature.name]] ) == TRUE, ] ) );
    
  # Remove genes which 'ENTREZID' could not be find
  tmp$signature <- tmp$signature[ ! is.na(tmp$signature$entrezgene), ];

  # Get features weights. If signature has no expression information set weights to 1
  tmp$signature$feature.weights <- ifelse(
    test = signatures$data[ tmp$signature$hgnc_symbol, signature.name ] == "0",
    yes  = rep( 1, nrow(tmp$signature) ),
    no   = signatures$data[ tmp$signature$hgnc_symbol, signature.name ]
    );
  
  # Score metabric samples using signature
  tmp$metabric.signature.scores <- score.multifeature.signature(
    expression.data = metabric$all.data, 
    survival.data   = metabric$all.survobj, 
    feature.names   = paste0( tmp$signature$entrezgene, "_at" ),
    feature.weights = tmp$signature$feature.weights
    );
   
  # Normalize samples' scores between -1 and 1
  # Change samples dichotomization from median of scores to percentile of events
  tmp$metabric.signature.scores <- data.frame(
    normalized.score = tmp$metabric.signature.scores$scores / nrow(tmp$signature),
    sample.group     = tmp$metabric.signature.scores$scores > quantile( tmp$metabric.signature.scores$scores, metabric.analysis$percentile.of.events ),
    row.names        = as.character( tmp$metabric.signature.scores$patient.IDs )
    );
  
  metabric.analysis$sample.scores[[signature.name]] <- tmp$metabric.signature.scores[ rownames(metabric.analysis$sample.scores), "normalized.score" ];

  # Create Kaplan-Meier plots
  metabric.analysis$km.plots[[signature.name]] <- create.km.plot(
    main = gsub("_", " ", signature.name),
    survival.object = metabric.analysis$surv.obj,
    patient.groups = factor( tmp$metabric.signature.scores$sample.group, levels=c(FALSE, TRUE) ),
    statistical.method = 'cox',
    
    main.cex = 1.5,
    
    key.stats.corner = c(0,0),
    key.stats.x.pos = 0.01,
    key.stats.y.pos = 0.02,
    key.stats.cex = 0.85,
        
    key.groups.corner = c(1,1),
    key.groups.x.pos = 0.95,
    key.groups.y.pos = 0.95,
    key.groups.cex = 0.85,
    
    xaxis.cex = 0.75,
    xlab.label = 'Time (Years)',
    xlab.cex = 1,

    yaxis.cex = 0.75,    
    ylab.label = 'Estimated Survival Proportion',
    ylab.cex = 1,
    
    enable.warnings = FALSE,
    key.groups.labels = c("Low", "High"),
    show.risktable = FALSE,
    
    # Adds margins to the plots
    top.padding = 1,
    right.padding = 1,
    bottom.padding = 1,
    left.padding = 1
    );

  # Calculate prognostic value for signature
  tmp$sign.odds.ratios <- calcHazardRatios( 
    surv.obj = metabric.analysis$surv.obj, 
    sample.group = tmp$metabric.signature.scores$sample.group, 
    signature.name = signature.name
    );

  metabric.analysis$sign.odds.ratios <- rbind(
    metabric.analysis$sign.odds.ratios, 
    tmp$sign.odds.ratios[ colnames(metabric.analysis$sign.odds.ratios) ]
    );
  }

# metabric.analysis$sign.odds.ratios[ metabric.analysis$sign.odds.ratios$pval < 0.05, ];

# Calculate statistics for each sample
metabric.analysis$sample.metascores <- data.frame( 
  mean = apply( metabric.analysis$sample.scores, 1, mean ),
  median = apply( metabric.analysis$sample.scores, 1, median ),
  sd = apply( metabric.analysis$sample.scores, 1, sd )
  );

# Sort samples by median score of all angiogenesis signatures and resort other datasets
tmp$sample.order <- order( metabric.analysis$sample.metascores$median );
metabric.analysis$sample.metascores <- metabric.analysis$sample.metascores[ tmp$sample.order, ];
metabric.analysis$sample.scores <- metabric.analysis$sample.scores[ tmp$sample.order, ];
metabric.analysis$survival <- metabric.analysis$survival[ tmp$sample.order, ];
metabric.analysis$surv.obj <- metabric.analysis$surv.obj[ tmp$sample.order ];

# Add x-coordenates for plot
metabric.analysis$sample.metascores$x <- 1:nrow( metabric.analysis$sample.metascores );

# Create top plot as scatterplot (whiskers == SD)
metabric.analysis$metascores.plot <- create.scatterplot(
  #
  formula = median ~ x,
  data = metabric.analysis$sample.metascores,
  
  xlimits = c(1, nrow(metabric.analysis$sample.metascores) ),
  ylimits = c(-1,1),
  
  abline.h = 0,
  abline.col = "darkgray",
  abline.lwd = 0.5,
  abline.lty = 2,
  
  cex = 0.2,
  
  # Specifying error bars
  error.bar.lwd = 0.2,
  error.whisker.angle = 0,
  y.error.up = metabric.analysis$sample.metascores$sd,
  y.error.bar.col = alpha("darkgray", 0.5),
  
  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Define order of rows in heatmap (chronological)
tmp$study.order <- c( "Hu_2004", "Mendiola_2008", "Bentink_2012", 
                      "Masiero_2013", "Khong_2013", "Pinato_2013", 
                      "Sanmartin_2014", "Langlois_2014", "Stefansson_2015" );

# Create heatmap of scores for angiogenic signatures
metabric.analysis$signatures.heatmap <- create.heatmap(
  
  x = metabric.analysis$sample.scores[, rev(tmp$study.order) ],
  
  # Disable clustering of rows and columns. Axis should be displayed as are in data frame (samples median score x study year)
  cluster.dimensions = "none",
  
  # Use rownames and colnames as labels to xaxis and yaxis
  xaxis.lab = NULL,
  yaxis.lab = gsub("_", " ", rev(tmp$study.order) ),
  
  description = 'Heatmap created using BoutrosLab.plotting.general'
  );

# Creat a heatmap to be used as colour key for signature.heatmap
metabric.analysis$key.heatmap <- create.heatmap(
  x = data.frame( x = -50:50 ),
  clustering.method = 'none',
  scale.data = FALSE,
  print.colour.key = FALSE,
  yaxis.tck = 0
  );

# Create forest plot of hazard ratios for each angiogenesis signature
# Transform values in log2
metabric.analysis$sign.odds.ratios$log2.HR <- log2(metabric.analysis$sign.odds.ratios$HR);
metabric.analysis$sign.odds.ratios$log2.lower_95CI <- log2(metabric.analysis$sign.odds.ratios$lower_95CI);
metabric.analysis$sign.odds.ratios$log2.upper_95CI <- log2(metabric.analysis$sign.odds.ratios$upper_95CI);

# Calculate y-axis position at forestplot by matching sample order of signatures heatmap
metabric.analysis$sign.odds.ratios$y <- match( 
  rownames( metabric.analysis$sign.odds.ratios ),
  gsub(" ", "_", metabric.analysis$signatures.heatmap$y.scales$labels)
  );

# Create forest plot
metabric.analysis$forestplot <- create.scatterplot(
  formula = y ~ log2.HR,
  data = metabric.analysis$sign.odds.ratios,
  pch = 15,

  # Add vertical line at HR 1
  abline.v = 0,
  abline.col = "darkgray",
  abline.lty = 3,
  
  # Specifying error bars
  error.bar.length = 0,
  x.error.left = metabric.analysis$sign.odds.ratios$log2.HR - metabric.analysis$sign.odds.ratios$log2.lower_95CI,
  x.error.right = metabric.analysis$sign.odds.ratios$log2.upper_95CI - metabric.analysis$sign.odds.ratios$log2.HR,

  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Build one metasignature from consensus of the angiogenesis signatures
# barplot( table( rowSums( sign( metabric.analysis$sample.scores ) ) ) ) 
metabric.analysis$consensus.metasignature$data <- data.frame( 
  # consensus.metasignature = sign( rowSums( sign( metabric.analysis$sample.scores ) ) ),
  consensus.metasignature = dichotomize.dataset( rowSums( sign( metabric.analysis$sample.scores ) ) ),
  row.names = names( rowSums( sign( metabric.analysis$sample.scores ) ) )
  );

# Create heatmap 
metabric.analysis$consensus.metasignature$heatmap <- create.heatmap(
  x = metabric.analysis$consensus.metasignature$data,
  clustering.method = 'none',
  
  # Use rownames and colnames as labels to xaxis and yaxis
  yat = 1,
  yaxis.lab = "Consensus",
  
  print.colour.key = FALSE
  );

# >>>>>>>>> MEDIAM DICHOTOMIZATION DOES NOT WORK PROPERLY WHEN MEDIAM HAS TIES <<<<<<<<<<<<<<<<<<<<<<<
# table( metabric.analysis$consensus.metasignature$data ) 
# dichotomize.dataset( c(-1,-1,-1, 0, 0, 0, 1, 1) )

# Create forest plot of hazard ratios for consensus signature
# Calculate prognostic value of consensus metasignature
metabric.analysis$consensus.metasignature$odds.ratios <- calcHazardRatios( 
  surv.obj = metabric.analysis$surv.obj, 
  sample.group = metabric.analysis$consensus.metasignature$data$consensus.metasignature, 
  signature.name = "Consensus"
  );

# Transform values in log2
metabric.analysis$consensus.metasignature$odds.ratios$log2.HR <- log2(metabric.analysis$consensus.metasignature$odds.ratios$HR);
metabric.analysis$consensus.metasignature$odds.ratios$log2.lower_95CI <- log2(metabric.analysis$consensus.metasignature$odds.ratios$lower_95CI);
metabric.analysis$consensus.metasignature$odds.ratios$log2.upper_95CI <- log2(metabric.analysis$consensus.metasignature$odds.ratios$upper_95CI);

# Create forest plot
metabric.analysis$consensus.metasignature$forestplot <- create.scatterplot(
  formula = 1 ~ log2.HR,
  data = metabric.analysis$consensus.metasignature$odds.ratios,

  # Use same xlimits of signatures' forestplot (avoid error during plot)
  xlimits = metabric.analysis$forestplot$x.limits, 
  ylimits = c(0.5,1.5),
  pch = 15,

  # Add vertical line at HR 1
  abline.v = 0,
  abline.col = "darkgray",
  abline.lty = 3,
  
  # Specifying error bars
  error.bar.length = 0,
  x.error.left = metabric.analysis$consensus.metasignature$odds.ratios$log2.HR - metabric.analysis$consensus.metasignature$odds.ratios$log2.lower_95CI,
  x.error.right = metabric.analysis$consensus.metasignature$odds.ratios$log2.upper_95CI - metabric.analysis$consensus.metasignature$odds.ratios$log2.HR,

  description = 'Scatter plot created by BoutrosLab.plotting.general'
  );

# Define breaks used in forest plots
tmp$forest.breaks <- c( 0.66, 1, 1.5 );

### Create multiplot with all plots
metabric.analysis$multiplot <- create.multiplot(
  #
  plot.objects = list( metabric.analysis$key.heatmap, 
                       metabric.analysis$consensus.metasignature$heatmap, metabric.analysis$consensus.metasignature$forestplot, 
                       metabric.analysis$signatures.heatmap, metabric.analysis$forestplot, 
                       metabric.analysis$metascores.plot ),
  
  panel.heights = c(0.3, 1, 0.08 ,0.05),
  panel.widths = c(1, 0.1),
  plot.layout = c(2, 4),
  layout.skip = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),

  axes.lwd = 0.5,
  
  xaxis.fontface = 1,
  xaxis.cex = 0.5,
  xaxis.tck = c(0.25, 0.25),
  xlab.label = c('\t', '\t', '        Samples', '\t', '\t', '            HR [95% CI]'),
  xlab.cex = 0.7,
  xlab.to.xaxis.padding = -1,
  x.spacing = -0.75,
  xlimits = list(
    NULL,
    NULL,
    0.8 * c(-1,1),
    NULL,
    0.8 * c(-1,1),
    NULL
    ),
  xat = list(
    c(1,51,101),
    NULL,
    log2( tmp$forest.breaks ),
    NULL,
    log2( tmp$forest.breaks ),
    NULL
    ),
  xaxis.lab = list(
    c(-1,0,1),
    NULL,
    tmp$forest.breaks,
    NULL,
    NULL,
    NULL
    ),

  yaxis.fontface = 1,
  yaxis.cex = 0.5,
  yaxis.tck = c(0.25, 0.25),
  ylab.label = c('Median score' , '\t', 'Angiogenesis signatures', '\t', '\t'),
  ylab.cex = 0.7,
  ylab.padding = 7,
  y.spacing = c(-1,-1,-1),
  ylimits = list(
    c(0.6,1),
    c(0.6,1),
    NULL,
    NULL,
    c(0, nrow(metabric.analysis$sign.odds.ratios) )+0.5,
    c(-1,1)
    ),
  yat = list(
    NULL,
    mean( c(0.6,1) ),
    NULL,
    1:length(metabric.analysis$signatures.heatmap$y.scales$labels),
    NULL,
    seq(-1,1,0.5)
    ),
  yaxis.lab = list(
    NULL,
    metabric.analysis$consensus.metasignature$heatmap$y.scales$labels,
    NULL,
    metabric.analysis$signatures.heatmap$y.scales$labels,
    NULL,
    seq(-1,1,0.5)
    ),

  print.new.legend = TRUE
  );

# Save plot on a file
tiff(
  filename = "aux_files/plots/angio_signatures_on_metabric.tiff", 
  type = "cairo", 
  width = 17, 
  height = 10, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

plot( metabric.analysis$multiplot );

dev.off();

#### Save KM-plots on a file
tiff(
  filename = "aux_files/plots/km_curves.tiff", 
  type = "cairo", 
  width = 30, 
  height = 25, 
  units = 'cm', 
  res = 500, 
  compression = 'lzw'
  );

# Create merged images
grid.arrange( 
  grobs = metabric.analysis$km.plots[tmp$study.order],
  layout_matrix = matrix( 
    data = 1:length(metabric.analysis$km.plots),
    nrow = length(metabric.analysis$km.plots)/3,
    ncol = 3,
    byrow = TRUE
    )
  );

dev.off();

### WRITE SESSION PROFILE TO FILE #################################################################
save.session.profile( 'aux_files/signatures_comparison_session_info.txt' );

save.image( file = "aux_files/signatures_comparison.Rdata" );
