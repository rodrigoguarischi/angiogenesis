#! /usr/bin/env Rscript

### null_model.R ##################################################################################
# This script was used to create the null distribution of rf models using Metabric dataset

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd( "~/angiogenesis_signature/" );

# Load required libraries and datasets
library( BoutrosLab.statistics.survival );
library( randomForestSRC );

# ### Only for construction of null.distribution object 
# 
# library( org.Hs.eg.db );
# library( BoutrosLab.datasets.breast.cancer );
# 
# ### Function - factor2numeric() #######################################################################
# # Input variables: A factor with values to be converted
# #
# # Output variables: A vector with numeric values
# #
# # Description: Convert a factor to a vector of numeric values
# factor2numeric <- function( x ) {
#   
#   stopifnot( inherits( x, "factor" ) );
#   
#   return( as.numeric( levels(x) )[x] );
#   }

### DATA ANALYSIS #################################################################################

# ### Code used to prepare the null.distribution object for later use
# 
# # Load Metabric dataset
# metabric <- BoutrosLab.datasets.breast.cancer::load.breast.cancer.datasets( datasets.to.load=c('Metabric.Training','Metabric.Validation') );
# 
# # Initialize object 'null.distribution' that will hold all data
# null.distribution <- list(); 
# 
# # Initialize object 'tmp' that will hold all temporary data
# tmp <- list();
# 
# # Obtain expression data
# tmp$training.continuous <- data.frame( t( metabric$all.data$Metabric.Training ) );
# tmp$validation.continuous <- data.frame( t( metabric$all.data$Metabric.Validation ) );
# 
# null.distribution$entrezIDs$training <- gsub("X", "", gsub("_at", "", colnames(tmp$training.continuous) ) );
# null.distribution$entrezIDs$validation <- gsub("X", "", gsub("_at", "", colnames(tmp$training.continuous) ) );
# 
# # Replace EntrezIDs for GeneSymbols on downstream analysis
# colnames(tmp$training.continuous) <- mapIds(
#   x = org.Hs.eg.db,
#   keys = null.distribution$entrezIDs$training,
#   keytype="ENTREZID",
#   column = "SYMBOL",
#   multiVals="asNA"
#   );
# 
# colnames(tmp$validation.continuous) <- mapIds(
#   x = org.Hs.eg.db,
#   keys = null.distribution$entrezIDs$validation,
#   keytype="ENTREZID",
#   column = "SYMBOL",
#   multiVals="asNA"
#   );
# 
# # Use safe characters as colnames and avoid downstream problems
# colnames(tmp$training.continuous) <- make.names( colnames(tmp$training.continuous) );
# colnames(tmp$validation.continuous) <- make.names( colnames(tmp$validation.continuous) );
# 
# # Convert variables to binary by median-dichotomizing expression values (Vinnie's paper)
# tmp$training.binary <- sapply( tmp$training.continuous, function(x) factor( as.numeric( x > median(x) ) ) );
# tmp$validation.binary <- sapply( tmp$validation.continuous, function(x) factor( as.numeric( x > median(x) ) ) );
# 
# # Append '__B' to column names
# colnames(tmp$training.binary) <- paste0( colnames( tmp$training.binary ), "__B" );
# colnames(tmp$validation.binary) <- paste0( colnames( tmp$validation.binary), "__B" );
# 
# # Combine continuous and binary information
# null.distribution$training <- cbind( tmp$training.continuous, tmp$training.binary );
# null.distribution$validation <- cbind( tmp$validation.continuous, tmp$validation.binary );
# 
# # Add clinical variables "age"
# null.distribution$training$age <- metabric$all.annotation$Metabric.Training[ rownames( null.distribution$training ), "age" ];
# null.distribution$validation$age <- metabric$all.annotation$Metabric.Validation[ rownames( null.distribution$validation), "age" ];
# 
# # Add clinical variables "stage" (Transform nulls into NAs. Will throw some warnings, that's ok!)
# null.distribution$training$stage <- factor( factor2numeric( metabric$all.annotation$Metabric.Training[ rownames( null.distribution$training ), "stage" ] ) );
# null.distribution$validation$stage <- factor( factor2numeric( metabric$all.annotation$Metabric.Validation[ rownames( null.distribution$validation ), "stage" ] ) );
# 
# # Truncate survival at 15 years follow-up time
# null.distribution$surv.obj$training <- survival.truncate( list( metabric$all.survobj$Metabric.Training ), 15)[[1]];
# null.distribution$surv.obj$validation <- survival.truncate( list( metabric$all.survobj$Metabric.Validation ), 15)[[1]];
# 
# # Adds survival info
# null.distribution$training$time <- null.distribution$surv.obj$training[,1];
# null.distribution$validation$time <- null.distribution$surv.obj$validation[,1];
# null.distribution$training$status <- null.distribution$surv.obj$training[,2];
# null.distribution$validation$status <- null.distribution$surv.obj$validation[,2];
# 
# # Impute missing data
# null.distribution$training <- impute.rfsrc(
#   formula = Surv(time, status) ~ .,
#   data = null.distribution$training
#   );
# 
# null.distribution$validation <- impute.rfsrc(
#   formula = Surv(time, status) ~ .,
#   data = null.distribution$validation
#   );
# 
# # Save object for later use
# save( null.distribution, file = "null_distribution/null_distribution.RData");

### Code used to create and test random signatures

# Load null.distribution object
load( file = "./null_distribution/null_distribution.RData" );

# Get all angiogenesis features
null.distribution$angiogenesis.features <- as.character( read.table( "null_distribution/angiogenesis_features.txt" )[, "feature.name"] );

# Data frame data will store results
null.distribution$all.results <- data.frame(
  error.rate = numeric(), 
  p.value = numeric(),
  features = character()
  );

# Read command-line parameters
args <- commandArgs(TRUE);

if( length(args) < 2 ){
  stop( "Missing argument. Rscript null_model.R <signature_size> <full|subset> " );
  }

# Load function that helps with assigning classes to samples
source( "predict_class.R" );

# Build 10k random forest models (null distribution)
for ( i in 1:10000 ) {

  # Define signature size
  null.distribution$signature.size <- as.numeric( args[1] );
  
  if ( tolower( args[2] ) == "full" ) {
    null.distribution$selected.features <- sample( 
      x = colnames( null.distribution$training )[! colnames(null.distribution$training) %in% c("time", "status") ],
      size = null.distribution$signature.size,
      replace = FALSE
      );
    } else if ( tolower( args[2] ) == "subset" ) {
    null.distribution$selected.features <- sample(
      x = null.distribution$angiogenesis.features[ ! null.distribution$angiogenesis.features %in% c("time", "status") ],
      size = null.distribution$signature.size,
      replace = FALSE
      );
    } else {
      stop( paste0("Invalid argument '", args[2] , "', must be full (all features on array) or subset (angiogenesis only features)" ) );
    }

  # Build a model
  null.distribution$fit <- rfsrc(
    formula = Surv(time, status) ~ .,
    data = null.distribution$training[, c(null.distribution$selected.features, "time", "status") ],
    
    # All data should already have been imputed but, just in case...
    na.action = "na.impute",
    
    # Options to speed up processing time
    importance = FALSE,
    tree.err = FALSE,
    nsplit = 10,
    ntime = 10
    );

  # Run prediction
  null.distribution$prediction <- predict.rfsrc( 
    object = null.distribution$fit, 
    newdata = null.distribution$validation,
    
    # All data should already have been imputed but, just in case...
    na.action = "na.impute",
    
    # Options to speed up processing time
    importance = FALSE
    );
  
  # Assign a class for each sample based on the probability that the patient has survived
  # low-risk          = p(survival) > 0.5 at 15 years
  # intermediate-risk = p(survival) > 0.5 at 7.5 years
  # high-risk         = p(survival) < 0.5 at 7.5 years
  null.distribution$predicted.classes <- predict.class(
      rfsrc.obj = null.distribution$prediction, 
      breaks = c(7.5, 15)
      );

  # Set p.value to 1 if all samples were assigned to only one class. logrank test would fail
  if( all( summary( null.distribution$predicted.classes ) != nrow(null.distribution$validation) ) ){

    # Calculate p-value for the model using logrank test
    null.distribution$p.value <- logrank.analysis(
      survival.object = null.distribution$surv.obj$validation,
      groups = null.distribution$predicted.classes
     )$pvalue[1];
    
    } else {
      
    null.distribution$p.value <- 1;
    
    }

  null.distribution$model.result <- data.frame(
    error.rate = null.distribution$prediction$err.rate[ length(null.distribution$prediction$err.rate) ],
    p.value = null.distribution$p.value,
    features = paste(null.distribution$selected.features, collapse = ",")
    );
  
  # Append model's result to data.frame
  null.distribution$all.results <- rbind(null.distribution$all.results, null.distribution$model.result );
  
  }

# Assign a filename to the output following the pattern <subset|full>_<signature_size>genes.csv
out.filename <- paste0(
  "null_distribution/results/",
  ifelse( tolower( args[2] ) == "full", "full_", "subset_" ),
  args[1], 
  "features.csv"
  );

# Write output file
write.csv(
  x = null.distribution$all.results,
  file = out.filename,
  row.names = FALSE
  );