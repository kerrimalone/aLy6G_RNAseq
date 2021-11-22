#################################
# Install the required packages # ----
#############################
libs_load <- function(x){
  for( i in x ){
    print(paste0("Checking for library: ", i))
    if(require( i , character.only = TRUE ) ){
      print(paste0(i, " already installed. Loading now"))
    }
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      print(paste0(i, " not installed. Trying CRAN for install."))
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
      paste0(i, " installed and loaded successfully")
    }
    if ( ! require(i, character.only=TRUE) ) {
      paste0(i," could not be installed from CRAN. Trying Bionconductor....")
      BiocManager::install(i)
      require( i , character.only = TRUE )
      paste0(i, " installed and loaded successfully")
    }
    if ( ! require(i, character.only=TRUE) ) {
      paste0(i, "could not be installed. Check manually")
    }
    #  Load package after installing
  }
}

libs_load(c("AnnotationFuncs","Biobase", "edgeR", "grDevices", 
            "grid", "gridExtra", "org.Mm.eg.db", "RColorBrewer", 
            "reshape2", "rtracklayer", "sigora", "slam", 
            "tidyverse", "tximport"))

`%notin%` <- Negate(`%in%`)
##########################
# Load required packages # ----
##########################



#######################
# Load data
#######################

# Summarise transcript abundance output from Salmon at the gene level using tximport
INPUT_FILE_DIR = "data/input/placenta"
GTF_FILE_OUT = paste(INPUT_FILE_DIR, "gencode.vM23.annotation.gtf.gz", sep="/")
GTF_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz"
GTF_FILE = paste(INPUT_FILE_DIR, "gencode.vM23.annotation.gtf", sep="/")
OUTPUT = "data/output/placenta"
PLOT_DIR = paste(OUTPUT, "plots", sep="/")
METADATA_FILE = paste0(INPUT_FILE_DIR, "/placenta_metadata.csv")


salmon_files <- list.files(INPUT_FILE_DIR, pattern = "quant", full.names = TRUE)
sample_names <- basename(salmon_files)
sample_names <- gsub("_quant", "", sample_names)

#Make annotation file
if(!file.exists(GTF_FILE)){
  download.file(GTF_URL, GTF_FILE)
  R.utils:::gunzip(GTF_FILE_OUT)
}


gtf <- rtracklayer::import(GTF_FILE)
gtf<-as.data.frame(gtf)
coding_tx = gtf %>% filter(type == "transcript") %>% filter(gene_type == "protein_coding")
tx_gene_df<- coding_tx %>% dplyr::select(transcript_id, gene_id)


# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
# (ignoreTxVersion = TRUE otherwise the cow transcripts names in the UCSC TxDb 
# won't match the ones in salmon's quant.sf)
mouse_txi <- tximport(salmon_files,
                   type = "salmon",
                   tx2gene = tx_gene_df)

colnames(mouse_txi$counts) <- sample_names
colnames(mouse_txi$abundance) <- sample_names



# Counts
mouse_txi$counts<-as.data.frame(mouse_txi$counts)
names(mouse_txi$counts)
anyDuplicated(rownames(mouse_txi$counts))
mouse_txi_counts<-(mouse_txi$counts)
head(mouse_txi_counts)
rownames(mouse_txi_counts)
mouse_txi_counts$gene_name<-rownames(mouse_txi_counts)

#Take each column and save it to a new text file named after the column with gene rownames kept.
COUNTS_DIR_OUT = paste(OUTPUT, "/counts", sep="")
dir.create(COUNTS_DIR_OUT)

for (i in seq(1,ncol(mouse_txi_counts))){
  if (i < 7) {
    count_df <-data.frame(mouse_txi_counts[,i])
    colnames(count_df) <- "counts"
    count_df$gene <-rownames(mouse_txi_counts)
    file_name <- paste(COUNTS_DIR_OUT, "/", colnames(mouse_txi_counts[i]), "_counts.txt", sep="")
    write.table(count_df, 
                file_name,
                append = FALSE,
                quote = FALSE,
                sep = "\t",
                row.names = TRUE,
                col.names = TRUE) 
  }
}
    


# Abundance/TPM
mouse_txi$abundance<-as.data.frame(mouse_txi$abundance)
names(mouse_txi$abundance)
anyDuplicated(rownames(mouse_txi$abundance))
mouse_txi_abundance<-(mouse_txi$abundance)
head(mouse_txi_abundance)
rownames(mouse_txi_abundance)
mouse_txi_abundance$gene_name<-rownames(mouse_txi_abundance)



#Take each column and save it to a new text file named after the column with gene rownames kept.
abundance_DIR_OUT = paste(OUTPUT, "/abundance", sep="")
dir.create(abundance_DIR_OUT)

for (i in seq(1,ncol(mouse_txi_abundance))){
  if (i < 7) {
    abundance_df <-data.frame(mouse_txi_abundance[,i])
    rownames(abundance_df) <-rownames(mouse_txi_abundance)
    colnames(abundance_df) <- "abundance"
    abundance_df$gene <- rownames(mouse_txi$abundance)
    file_name <- paste(abundance_DIR_OUT, "/", colnames(mouse_txi_abundance[i]), "_abundance.txt", sep="")
    write.table(abundance_df, 
                file_name,
                append = FALSE,
                quote = FALSE,
                sep = "\t",
                row.names = TRUE,
                col.names = TRUE) 
  }
}


count_files <- list.files(path = COUNTS_DIR_OUT,
                        pattern         = "*1_counts.txt",
                        all.files       = TRUE,
                        full.names      = FALSE,
                        recursive       = FALSE,
                        ignore.case     = FALSE)

# Reads and merges a set of files containing counts
counts <- readDGE(files = count_files, 
                  path= COUNTS_DIR_OUT,
                  header = TRUE, 
                  columns = c(2,1), 
                  sep="\t", 
                  comment.char = "#")

names(counts)
head(counts$samples)
head(counts$counts)
dim(counts)

# Output data
SAMPLES_FILE_OUT = paste(OUTPUT, "placenta_samples.txt", sep="/")
COUNT_FILE_OUT = paste(COUNTS_DIR_OUT, "placenta_counts.txt", sep="/")


write.table(x = counts$samples,file = SAMPLES_FILE_OUT, sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)
write.table(x = counts$counts,file = COUNT_FILE_OUT, sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)

#####################
# Clean input data
#####################

#Clean input files
raw.counts <- read.table(file = COUNT_FILE_OUT, header = TRUE)
colnames(raw.counts)
head(raw.counts)
#Remove extraneous file extenstion info 
colnames(raw.counts) <- gsub("_counts","",colnames(raw.counts))

head(raw.counts)

COUNT_FILE_OUT_CLEAN = paste(COUNTS_DIR_OUT, "placenta_samples_clean.txt", sep="/")

write.table(x = raw.counts,file = COUNT_FILE_OUT_CLEAN, sep = "\t",quote = FALSE, row.names = TRUE,col.names = NA)




# Create annotation table with counts information
annotated.counts<-raw.counts

#Remove versioning
rownames(annotated.counts)<-gsub("(.*)\\..*", "\\1", rownames(annotated.counts))


columns(org.Mm.eg.db)
as.list(org.Mm.egENSEMBL) %>% head()

# Get gene symbols from NCBI gene identifiers
annotated.counts$gene_name <- mapIds(org.Mm.eg.db,
                                     keys      = rownames(annotated.counts),
                                     column    = "GENENAME",
                                     keytype   = "ENSEMBL",
                                     multiVals = "first")

# Get ENSEMBL gene ids from NCBI gene identifiers
annotated.counts$ENSEMBL.tag <- mapIds(org.Mm.eg.db,
                                       keys      = rownames(annotated.counts),
                                       column    = "SYMBOL",
                                       keytype   = "ENSEMBL",
                                       multiVals = "first")

head(annotated.counts)
dim(annotated.counts)


COUNTS_FILE_OUT_CLEAN = paste(OUTPUT, "placenta_counts_clean_annotated.txt", sep="/")

# Output data
write.table(x         = annotated.counts,
            file      = COUNTS_FILE_OUT_CLEAN,
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)



##############
# DGElist
##############
# Create DGElist containing information about condition and annotation
head(raw.counts)

metadata<-read.table(file=METADATA_FILE,
                     sep = ",",header=TRUE)
head(metadata)
colnames(metadata)[1]<-"Sample"

condition <- relevel(factor(metadata$Group), ref = "iso_preg_placenta")
head(condition)

names(annotated.counts)
annotated.counts$gene_symbol<-rownames(annotated.counts)
gene.annotation <- dplyr::select(annotated.counts,
                                 gene_name,
                                 gene_symbol,
                                 ENSEMBL.tag)

head(gene.annotation)
dim(gene.annotation)
rownames(gene.annotation)

dim(raw.counts)
dim(gene.annotation)
head(raw.counts)
raw.counts$gene_symbol<-rownames(raw.counts)
raw.counts_clean<-merge(gene.annotation,raw.counts, by="gene_symbol")

#Should equal gene annotation that has removed duplicates etc.
dim(annotated.counts)
annotated.counts.clean<-annotated.counts[,1:6]
head(annotated.counts.clean)

alyg6_dgelist <- DGEList(counts       = annotated.counts.clean,
                         group        = condition,
                         genes        = gene.annotation,
                         lib.size     = NULL,
                         norm.factors = NULL,
                         remove.zeros = FALSE)

names(alyg6_dgelist)
dim(alyg6_dgelist)
head(alyg6_dgelist$counts)
dim(alyg6_dgelist$samples)
head(alyg6_dgelist$genes)
rownames(alyg6_dgelist$samples)



# Add metadata to DGElist$samples
head(metadata)
dim(metadata)
rownames(metadata)<-metadata[,1]
metadata<-metadata[,2:6]
colnames(metadata[1]) = "group"
rownames(metadata)<-gsub("NG\\-", "NG\\.", rownames(metadata))
alyg6_dgelist$samples <- merge(x  = alyg6_dgelist$samples,
                               y  = metadata,
                               by = "row.names") 

head(alyg6_dgelist$samples)
rownames(alyg6_dgelist$samples)<-alyg6_dgelist$samples[,1]
alyg6_dgelist$samples<-alyg6_dgelist$samples[,2:9]
names(alyg6_dgelist$samples)[1]<-"Sample"



##############################################
# Filtering of zero and lowly expressed tags # ----
##############################################
# Filter non expressed tags (all genes that have zero counts in all samples
#21756 genes
dim(alyg6_dgelist)
alyg6_no_zeros <-  alyg6_dgelist[!rowSums(cpm(alyg6_dgelist$counts) < 1) >=3, ]

dim(alyg6_no_zeros$counts)

head(alyg6_no_zeros$counts)



########################################################
# Normalization of data using Trimmed Mean of M-values # ----
#    (based on RNA composition between libraries)      #
########################################################

# Calculate normalisation factor for our DGElist.
# With edgeR, counts are not transformed in any way after normalization,
# instead normalization will modify library size.
alyg6_norm <- calcNormFactors(alyg6_no_zeros, method = "TMM")
head(alyg6_norm$samples) #check norm.factors before and after normalisation
# save.image("alyg6.Rdata")
#####################################################################
# Quality check of filtered libraries by plotting density of counts # ----
#####################################################################

DENSITY_FILT_PLOT = paste(PLOT_DIR, "Density_filt_placenta.png", sep="/")
# Log10 transform the filtered count data for better visualization
count_filt_log10 <- log10(alyg6_norm$counts[, 1 : ncol(alyg6_norm$counts)] + 1)

# Plot density of count for all libraries
png(filename = DENSITY_FILT_PLOT,
    width    = 1366,
    height   = 768,
    units    = "px")

plot(density(count_filt_log10[, 1]),
     main = "Density plot of count per gene post filtering",
     lty  = 1,
     xlab = "Log10 of count per gene",
     ylab = "Density",
     col  = "black",
     ylim = c(0.0,0.6))

for (i in 2 : ncol(count_filt_log10)) {
  lines(density(count_filt_log10[, i]),
        lty = 1,
        col = "black")
}

dev.off()

#############################################################
# Exploratory data analysis: Multidimensional scaling plots # ----
#############################################################
MDS_PLOT = paste(PLOT_DIR, "placenta_MDS.png", sep="/")

# Plot MDS of all samples
png(filename = MDS_PLOT,
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(alyg6_norm, labels=alyg6_norm$samples[,1])

dev.off()

# Plot MDS of all samples by animal
MDS_ANIMAL_PLOT = paste(PLOT_DIR, "placenta_MDS_animal.png", sep="/")

png(filename = MDS_ANIMAL_PLOT,
    width    = 1366,
    height   = 768,
    units    = "px")

plotMDS(alyg6_norm, labels=alyg6_norm$samples[,8])

dev.off()


placenta<-alyg6_norm

placenta_groups<-factor(placenta$samples$group)
placenta_groups<-as.factor(placenta_groups)
placenta_matrix <- model.matrix(~Condition.2,
                                  data = placenta$samples)
colnames(placenta_matrix)
dim(placenta_matrix)


DE_ANALYSIS_OUT = paste(OUTPUT, "DE_analysis", sep="/")
dir.create(DE_ANALYSIS_OUT)
placenta_MATRIX = paste(DE_ANALYSIS_OUT, "placenta_test_matrix.txt", sep ="/")
write.table(x         = placenta_matrix,
            file      = placenta_MATRIX,
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

############################################################################
# Estimate the dispersion parameter for each tag using the Cox-Reid method # ----
#                      (for multi-factor data)                             #
############################################################################
placentadisp <- estimateGLMCommonDisp(y    = placenta,
                                        design  = placenta_matrix,
                                        verbose = TRUE)

placentadisp <- estimateGLMTrendedDisp(y   = placentadisp,
                                         design = placenta_matrix)

placentadisp <- estimateGLMTagwiseDisp(y   = placentadisp,
                                         design = placenta_matrix)

names(placentadisp)


# Plot the dispersion

placenta_DISP_PLOT = paste(PLOT_DIR, "placenta_disp.png", sep="/")

png(filename = placenta_DISP_PLOT ,
    width    = 1366,
    height   = 768,
    units    = "px")

plotBCV(placentadisp )

dev.off()

# Show the calculated dispersion
placentadisp$common.dispersion

# And show its square root, the coefficient of biological variation #http://seqanswers.com/forums/showthread.php?t=5591
#0.24265729
sqrt(placentadisp$common.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information 
Tagwisedisp_placenta <- cbind(placentadisp$genes, placentadisp$tagwise.dispersion)
dim(Tagwisedisp_placenta)

# Write into a table the calculated tagwise dispersion
write.matrix(x    = Tagwisedisp_placenta,
             file = "placenta_Tagwise_dispersion_test_matrix.txt",
             sep  = "\t")

##################################################################
# Determine differential expression using quasi-negative binomial GLMs # ----
##################################################################
# Fit a quasi-negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
placenta_fit <- glmFit(y = placentadisp, design = placenta_matrix)
names(placenta_fit)
colnames(placenta_fit$design)
dim(placenta_fit)


#################################
#     preg_placenta_qlf   #
################################
placenta_lrt <- glmLRT(placenta_fit, coef = 2)
placenta_tags <- topTags(object        = placenta_lrt,
                           n             = "inf",
                           adjust.method = "BH")
nrow(placenta_tags)

ALL_GENES_LOGFC = paste0(OUTPUT, "/placenta_all_genes_logFC.tsv")

#Results for all genes
write.table(x         = placenta_tags$table,
            file      = ALL_GENES_LOGFC,
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)

FDR_0.05_DE_placenta<- subset(placenta_tags$table, FDR < 0.05)

#357
nrow(FDR_0.05_DE_placenta)


PLACENTA_FDR_0.05 = paste0(OUTPUT, "/placenta_DE_FDR_0.05.tsv")

write.table(x         = FDR_0.05_DE_placenta,
            file      = PLACENTA_FDR_0.05,
            sep       = "\t",
            quote     = FALSE,
            row.names = TRUE,
            col.names = NA)




############################
# Pathway analysis: SIGORA #
############################
# require(devtools)
#install_version("sigora", version = "3.0.5", repos = "http://cran.us.r-project.org")

# install.packages("slam")
library(sigora)
library(slam)


SIGORA_OUT = paste(OUTPUT, "sigora", sep="/")
dir.create(SIGORA_OUT)


data(reaM)
get_sigora_pathways <-function(DE_GENE_LIST, GENES_COL_INDEX) {
  sig.data.genes<-(DE_GENE_LIST[,GENES_COL_INDEX])
  sig.data.genes<-na.omit(sig.data.genes)
  OUTFILE_PATH = paste(SIGORA_OUT, deparse(substitute(DE_GENE_LIST)), sep="/")
  OUTFILE = paste(OUTFILE_PATH, "sigora.tsv", sep="_")
  sigRes<-sigora(queryList=sig.data.genes,GPSrepo=reaM,level=4, saveFile = OUTFILE)
}


get_sigora_pathways(FDR_0.05_DE_placenta,3)



