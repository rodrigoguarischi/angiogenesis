### signature_creation.R ##########################################################################
# This script was used to create the angiogenic signature using Metabric dataset

### PREAMBLE ######################################################################################
# Settings for environment

# Define working directory
setwd("~/angiogenesis_transcriptome/scripts/");

# Add some locally installed libraries to path
.libPaths( "/home/rodrigoguarischi/R/x86_64-pc-linux-gnu-library/3.3" );

# Load required libraries and datasets
library( BoutrosLab.plotting.general );
library( BoutrosLab.datasets.breast.cancer );
library( BoutrosLab.plotting.survival );
library( randomForestSRC );
library( parallel );                                # Enables parallel processing
library( org.Hs.eg.db );
# library( cpcgene.utilities );
library( gridExtra );
# library( caret );
# library( pROC );                                # Used to create ROC curves

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

# Use up to N cores on parallel processing. 
# Memory usage grows linearly with # of cores. No big improvement in processing time above 15 cores
options(mc.cores = ifelse( detectCores() > 15, 15, detectCores() - 1) );

# Load Metabric dataset
metabric <- load.breast.cancer.datasets( datasets.to.load=c('Metabric.Training','Metabric.Validation') );

# Initialize object 'metabric.analysis' that will hold data for signatures evaluation on Metabric dataset
metabric.analysis <- list(); 

# Angiogenesis DE genes from RNA-seq. NCBI geneIDs from human homologous
metabric.analysis$signature.genes <- unique( 
  read.delim(
    file = "aux_files/rnaseq_human_homologous.txt",
    sep = "\t",
    header = FALSE,
    comment.char = "#")[,5]
    );

# Assess which genes have expression data available
metabric.analysis$probes.names <- paste0(metabric.analysis$signature.genes, "_at");
metabric.analysis$mapped.genes <- metabric.analysis$probes.names[ metabric.analysis$probes.names %in% rownames( metabric$all.data$Metabric.Training ) ];
metabric.analysis$unmapped.genes <- metabric.analysis$probes.names[! metabric.analysis$probes.names %in% rownames( metabric$all.data$Metabric.Training ) ];

if( length(metabric.analysis$unmapped.genes) > 0 ){
  warning( "Genes ", paste(metabric.analysis$unmapped.genes, collapse = ", ") , " could not be mapped (", round( length( metabric.analysis$unmapped.genes ) / length( metabric.analysis$probes.names ) * 100 , 2), "% of total)." );
  }

# Subset expression data for mapped probes
metabric.analysis$training.continuous <- data.frame( t( metabric$all.data$Metabric.Training[ metabric.analysis$mapped.genes , ] ) );
metabric.analysis$validation.continuous <- data.frame( t( metabric$all.data$Metabric.Validation[ metabric.analysis$mapped.genes, ] ) );

# Replace EntrezIDs for GeneSymbols on downstream analysis
colnames(metabric.analysis$training.continuous) <- mapIds(
  x = org.Hs.eg.db,
  keys = gsub("X", "", gsub("_at", "", colnames(metabric.analysis$training.continuous) ) ),
  keytype="ENTREZID",
  column = "SYMBOL",
  multiVals="asNA"
  );
colnames(metabric.analysis$validation.continuous) <- mapIds(
  x = org.Hs.eg.db,
  keys = gsub("X", "", gsub("_at", "", colnames(metabric.analysis$validation.continuous) ) ),
  keytype="ENTREZID",
  column = "SYMBOL",
  multiVals="asNA"
  );

# Use safe characters as colnames and avoid downstream problems
colnames(metabric.analysis$training.continuous) <- make.names( colnames(metabric.analysis$training.continuous) );
colnames(metabric.analysis$validation.continuous) <- make.names( colnames(metabric.analysis$validation.continuous) );

# Convert variables to binary by median-dichotomizing expression values (Vinnie's paper)
metabric.analysis$training.binary <- sapply( metabric.analysis$training.continuous, function(x) factor( as.numeric( x > median(x) ) ) );
metabric.analysis$validation.binary <- sapply( metabric.analysis$validation.continuous, function(x) factor( as.numeric( x > median(x) ) ) );

# Append '__B' to column names
colnames(metabric.analysis$training.binary) <- paste0( colnames( metabric.analysis$training.binary ), "__B" );
colnames(metabric.analysis$validation.binary) <- paste0( colnames( metabric.analysis$validation.binary), "__B" );

# Combine continuous and binary information
metabric.analysis$training <- cbind( metabric.analysis$training.continuous, metabric.analysis$training.binary );
metabric.analysis$validation <- cbind( metabric.analysis$validation.continuous, metabric.analysis$validation.binary );

# Add clinical variables "age"
metabric.analysis$training$age <- metabric$all.annotation$Metabric.Training[ rownames( metabric.analysis$training ), "age" ];
metabric.analysis$validation$age <- metabric$all.annotation$Metabric.Validation[ rownames( metabric.analysis$validation), "age" ];

# Add clinical variables "stage" (Transform nulls into NAs. Will throw some warnings, that's ok!)
metabric.analysis$training$stage <- factor( factor2numeric( metabric$all.annotation$Metabric.Training[ rownames( metabric.analysis$training ), "stage" ] ) );
metabric.analysis$validation$stage <- factor( factor2numeric( metabric$all.annotation$Metabric.Validation[ rownames( metabric.analysis$validation ), "stage" ] ) );

# Truncate survival at 15 years follow-up time
metabric.analysis$surv.obj$training <- survival.truncate( list( metabric$all.survobj$Metabric.Training ), 15)[[1]];
metabric.analysis$surv.obj$validation <- survival.truncate( list( metabric$all.survobj$Metabric.Validation ), 15)[[1]];

# Adds survival info
metabric.analysis$training$time <- metabric.analysis$surv.obj$training[,1];
metabric.analysis$validation$time <- metabric.analysis$surv.obj$validation[,1];
metabric.analysis$training$status <- metabric.analysis$surv.obj$training[,2];
metabric.analysis$validation$status <- metabric.analysis$surv.obj$validation[,2];

# ### Save object to null distribution analysis
# 
# # Impute missing data
# metabric.analysis$training <- impute.rfsrc(
#   formula = Surv(time, status) ~ .,
#   data = metabric.analysis$training
#   );
# 
# metabric.analysis$validation <- impute.rfsrc(
#   formula = Surv(time, status) ~ .,
#   data = metabric.analysis$validation
#   );
# 
# # Save object 
# save( metabric.analysis, file = "null_distribution/metabric_analysis_obj.RData");

### Build signature

# Run univariate analysis to get significant features only
metabric.analysis$univariate.cox.results <- sapply(
  colnames(metabric.analysis$training)[ ! colnames(metabric.analysis$training) %in% c("time", "status") ],
  function (x) {
    fit.coxmodel(
      survobj = metabric.analysis$surv.obj$training,
      groups = metabric.analysis$training[[x]]
      )
    }
  );

metabric.analysis$univariate.cox.results <- data.frame( t(metabric.analysis$univariate.cox.results) );
colnames(metabric.analysis$univariate.cox.results) <- c("HR", "HR_min", "HR_max", "p", "n");

# Order features by p-values
metabric.analysis$univariate.cox.results <- metabric.analysis$univariate.cox.results[ order(metabric.analysis$univariate.cox.results$p), ];

# Filter only significant features ( p < 0.05 )
metabric.analysis$sig.features <- rownames( metabric.analysis$univariate.cox.results[ metabric.analysis$univariate.cox.results$p < 0.05, ] );

# "Stage" gives a infinite beta. Forcing inclusion of stage feature
metabric.analysis$sig.features <- c( "stage", metabric.analysis$sig.features ); 

### Model creation

# Create object 'model' to store model infos
model <- list();

# Incrementally try features to get model with minimum error
for ( n.var in 2:length(metabric.analysis$sig.features) ){

  print( paste0("Using top ", n.var , " features" ) );
  
  # Select top N features
  selected.features <- metabric.analysis$sig.features[1:n.var];
  
  model$fit[[paste0("top", n.var)]] <- rfsrc(
    formula = Surv(time, status) ~ .,
    data = metabric.analysis$training[, c(selected.features, "time", "status") ],

    # Impute missing data
    na.action = "na.impute",
    
    # Save forest metrics
    importance = TRUE,
    tree.err = TRUE
    );

  }

# Calculate OOB error for each model
model$oob.error <- sapply( model$fit, function(x) x$err.rate[ length( x$err.rate ) ] );

# Get name and used predictors of model with minimum error
model$best.model$name <- names( which( model$oob.error == min( model$oob.error ) ) );
model$fit$final <- model$fit[[model$best.model$name]];

# Plot model accucary versus mtrys
png(
  filename = "plots/model_accuracy.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = "white",
  res = 100
  );

# Plot model variables importance
png(
  filename = "~/Desktop/variables_importance.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 12,
  pointsize = 12,
  bg = "white",
  res = 100
  );
plot( varImp( model$fit$final, scale = FALSE ) );
dev.off();

# Plot density plots for features
tiff(
  filename = "aux_files/plots/features_density.tiff",
  type = "cairo",
  width = 30,
  height = 25,
  units = 'cm',
  res = 500,
  compression = 'lzw'
  );
# png(
#   filename = "aux_files/plots/features_density.png",
#   type = "cairo",
#   units = "cm",
#   width = 30,
#   height = 25,
#   bg = "white",
#   res = 300
#   );

plot.variable(model$fit$final, surv.type = "rel.freq", cex.lab = 1.6, cex.axis = 1.45); 

dev.off();

#########################################################################################################
######################################### APPLY MODEL ###################################################
#########################################################################################################

# Predict sample classes
metabric.analysis$prediction <- predict.rfsrc(
  object = model$fit$final,                               # Feed the model into prediction step
  newdata = metabric.analysis$validation,                  # Use validation set
  na.action = "na.impute"
  );

# Load function that helps with assigning classes to samples
source( "predict_class.R" );

# Assign a class for each sample based on the probability that the patient has survived
# low-risk          = p(survival) > 0.5 at 15 years
# intermediate-risk = p(survival) > 0.5 at 7.5 years
# high-risk         = p(survival) < 0.5 at 7.5 years
metabric.analysis$validation$predicted.class <- predict.class( rfsrc.obj = metabric.analysis$prediction, breaks = c(7.5, 15) );

# # Fit a cox model
# fit.coxmodel(
#   survobj = metabric.analysis$surv.obj$validation,
#   groups = factor( metabric.analysis$validation$predicted.class, levels = c("low", "intermediate", "high") ),
#   return.cox.model = TRUE
#   );

# Plot a KM curve
png(
  filename = "plots/km_curve.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = "white",
  res = 100
  );
create.km.plot(
  survival.object = metabric.analysis$surv.obj$validation,
  patient.groups = factor( metabric.analysis$validation$predicted.class, levels = c("low", "intermediate", "high") ),
  xlab.label = 'Time (Years)'
  );
dev.off();


png(
  filename = "~/Desktop/km.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = NA,
  res = 100
  );




a <- metabric.analysis$validation

factor(a$stage)
a <- a[ a$stage %in% 1:4, ];
create.km.plot(
  survival.object = Surv( a$time, a$status ),
  patient.groups = a$stage,
  xlab.label = 'Time (Years)'
);
dev.off();

# Plot ROC curve
metabric.analysis$prediction.prob <- predict(
  model$fit$final,                                 # Feed the model into prediction step
  newdata = metabric.analysis$validation,          # Use validation set
  na.action = "na.impute",
  type = "prob"
  );

png(
  filename = "plots/roc.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = "white",
  res = 100
  );

roc(
  response = metabric.analysis$validation$survived,  
  predictor = metabric.analysis$validation$prediction.prob$High,
  
  # Create plot
  plot = TRUE,

  # Adds AUC (w/ confidence intervals) info to the plot
  print.auc = TRUE,
  ci = TRUE,

  # Identity line
  identity.lty=2,
  identity.lwd=2
  );
dev.off();

### Subset classification by PAM50 subtypes

# Pull PAM50 subtypes information
# metabric.analysis$validation$subtypes <- metabric$all.subtype$Metabric.Validation;
metabric.analysis$validation$subtypes <- metabric$all.annotation$Metabric.Validation[rownames(metabric.analysis$validation), "subtype"];

for( subtype in c("Basal", "Her2", "LumA", "LumB", "Normal") ){
  
  metabric.analysis$subtype.sample <- which( metabric.analysis$validation$subtypes == subtype );
  
  metabric.analysis$plot[[ paste0( subtype, "_survival" )]] <- create.km.plot(
    survival.object = metabric.analysis$surv.obj$validation[metabric.analysis$subtype.sample],
    patient.groups = factor( metabric.analysis$validation$predicted.class[metabric.analysis$subtype.sample], levels = c("low", "intermediate", "high") ),
    xlab.label = 'Time (Years)',
    main = subtype
    );

  }

png(
  filename = "aux_files/plots/PAM50subtypes_survival.png",
  type = "cairo",
  units = "in",
  width = 15,
  height = 20,
  pointsize = 12,
  bg = "white",
  res = 100
  );

# Create merged images
grid.arrange( 
  grobs = metabric.analysis$plot,
  as.table = TRUE
  );

dev.off();

### Save model and variables for latter use
# load(file = "aux_files/signature_creation.RData");
save( list = ls(), file = "aux_files/signature_creation.RData");

#########################################################################################################
############################# APPLY MODEL (Prostate Cancer set) #########################################
#########################################################################################################

# Initialize object 'cpcgene' that will hold data for signatures evaluation on CPC-GENE dataset
cpcgene <- list();

# Load pre-processed expression data file
cpcgene$array.data <- read.table(
  file = "/.mounts/labs/cpcgene/private/Microarrays/mRNA/output/batch_correct_annotation/2016-02-05_cpcgene_mRNA_array_combat_expression_GRCh38.txt", 
  sep = '\t',
  header = TRUE,
  row.names = 1,
  as.is = TRUE
  );

# Assess which genes have expression data available
# Get all angiogenesis genes even though only a subset will be used on the model
cpcgene$probes.names <- gsub("_at", "", metabric.analysis$probes.names );
cpcgene$mapped.genes <- cpcgene$probes.names[ cpcgene$probes.names %in% rownames( cpcgene$array.data ) ];
cpcgene$unmapped.genes <- cpcgene$probes.names[! cpcgene$probes.names %in% rownames( cpcgene$array.data ) ];

if( length(cpcgene$unmapped.genes) > 0 ){
  warning( "Genes ", paste(cpcgene$unmapped.genes, collapse = ", ") , " could not be mapped (", round( length( cpcgene$unmapped.genes ) / length( cpcgene$probes.names ) * 100 , 2), "% of total)." );
  }

# Subset angiogenesis expression data
cpcgene$mRNA.continuous <- data.frame( t( cpcgene$array.data[ cpcgene$mapped.genes, grepl( "^CPCG", colnames(cpcgene$array.data) ) ] ) );

# Replace EntrezIDs for GeneSymbols on downstream analysis
colnames(cpcgene$mRNA.continuous) <- mapIds(
  x = org.Hs.eg.db,
  keys = gsub("X", "", colnames(cpcgene$mRNA.continuous) ),
  keytype="ENTREZID",
  column = "SYMBOL",
  multiVals="asNA"
  );

# Use safe characters as colnames and avoid downstream problems
colnames(cpcgene$mRNA.continuous) <- make.names( colnames(cpcgene$mRNA.continuous) );

# Convert variables to binary by median-dichotomizing expression values (Vinnie's paper)
cpcgene$mRNA.binary <- sapply( cpcgene$mRNA.continuous, function(x) factor( as.numeric( x > median(x) ) ) );

# Append '__B' to column names
colnames(cpcgene$mRNA.binary) <- paste0( colnames( cpcgene$mRNA.binary ), "__B" );

# Combine continuous and binary information
cpcgene$mRNA.full <- cbind( cpcgene$mRNA.continuous, cpcgene$mRNA.binary );

# Import sample annotation data
cpcgene$ann <- get.patient.annotation( x = rownames( cpcgene$mRNA.full ) );

# Adds age at treatment 
cpcgene$mRNA.full$age <- cpcgene$ann$age_at_treatment;

# Adds T-category 
cpcgene$mRNA.full$t.category <- cpcgene$ann$clinical_t;
cpcgene$mRNA.full$stage <- factor( gsub("^.*(\\d+).*$", "\\1", cpcgene$mRNA.full$t.category) );

# Adds survival object (Biochemical Recurrency)
cpcgene$mRNA.full$surval.obj <- Surv( cpcgene$ann$time_to_bcr, cpcgene$ann$bcr );

# Drop samples without survival information
cpcgene$mRNA.full <- cpcgene$mRNA.full[ ! is.na( cpcgene$mRNA.full$surval.obj ), ];

# Truncate survival object at 10 years
cpcgene$mRNA.full$surval.obj <- survival.truncate( list( cpcgene$mRNA.full$surval.obj ), 10)[[1]];

# Apply model on CPC-GENE data
cpcgene$prediction <- predict.rfsrc(
  model$fit$final,                       # Feed model into prediction step
  newdata = cpcgene$mRNA.full            # Use validation set
  );

cpcgene$prediction$predicted.class <- predict.class( rfsrc.obj = cpcgene$prediction, breaks = c(5, 10) );

# # Fit a cox model
# fit.coxmodel(
#   survobj = cpcgene$mRNA.full$surval.obj,
#   groups = factor( cpcgene$mRNA.full$prediction, levels=c("Low", "High") ),
#   return.cox.model = TRUE
#   );

# Plot a KM curve
png(
  filename = "plots/km_curve.cpcgene.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = "white",
  res = 100
  );

create.km.plot(
  survival.object = cpcgene$mRNA.full$surval.obj,
  patient.groups = cpcgene$prediction$predicted.class,
  xlab.label = 'Time (Years)'
  );

dev.off();

#########################################################################################################
########### Evaluate model's performance relative to the null distribution ##############################
#########################################################################################################

# Run null_model.R script to create random signatures
# 
# cd ~/angiogenesis_signature/; qsub -N NullDist -b y -t 1-1000 -m beas -M Rodrigo.Sousa@oicr.on.ca -cwd -e null_distribution/errors_out/ -o /dev/null "module use /oicr/local/boutroslab/Modules/modulefiles; module use /oicr/local/boutroslab/Modules/modulefiles-deb8; module load R-BL/2016-07-22; Rscript null_model.R"
# cd ~/angiogenesis_signature/; qsub -N NullDist -b y -t 1-1000 -cwd -e null_distribution/errors_out/ -o /dev/null "module use /oicr/local/boutroslab/Modules/modulefiles; module use /oicr/local/boutroslab/Modules/modulefiles-deb8; module load R-BL/2016-07-22; Rscript null_model.R"
# 

# Create object to hold all data
null.distribution <- NULL;

for( c_file in list.files(path = "null_distribution/results", pattern = "*.csv", full.names = TRUE ) ){
  print( paste0("Reading file: ", c_file) );
  null.distribution <- c(null.distribution, read.csv(file = c_file)$p.value );
  }

# Apply logrank test to obtain the p-value for the model
final.model.pvalue <- logrank.analysis(
    survival.object = metabric.analysis$surv.obj$validation,
    groups = factor( metabric.analysis$validation$prediction, levels = c("Low", "High") )
  )$pvalue[1]

# Plot a KM curve
png(
  filename = "plots/null_dist.png",
  type = "cairo",
  units = "in",
  width = 10,
  height = 10,
  pointsize = 12,
  bg = "white",
  res = 100
  );

plot( 
  x = density(-log10( null.distribution ) ),
  main = "Metabric Validation",
  xlab = "-log10(p-value)"
  );

arrows( 
  x0 = -log10( final.model.pvalue ),
  y0 = 0.04,
  y1 = 0.01,
  lwd = 3,
  col = "red"
  );

dev.off();

# Get number of significant models at p < 0.05
table( null.distribution < 0.05 )["TRUE"];
# 
# TRUE 
# 9041255 
# 

# Get number of models with smaller p-values and respective proportion
table( null.distribution < final.model.pvalue )["TRUE"]; 
prop.table( table( null.distribution < final.model.pvalue ) )["TRUE"]*100;
# 
# TRUE 
# 49666 
# 
# TRUE 
# 0.49666 
# 