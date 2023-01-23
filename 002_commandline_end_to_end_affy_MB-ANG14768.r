#!/usr/bin/env Rscript
## this Rscript has command line arguments, it will check the sdfr file and select the appropriate "Characteristics" column to compare samples in limma. Filename and patient are always read.
 ####  Arguments: sdrf file; Column name to consider for DE; boolean for paired samples (optional)
#args = commandArgs(trailingOnly=TRUE)
print('Arguments: raw_data directory; sdrf file; Column name to consider for DiffExp; boolean for paired samples (optional)')
args = c("MB-ANG14768", 'MB-ANG14768.sdrf', 'Characteristics[group]', "False")

print(args)
print('this program is supposed to be run after python script 001')

############# Installed packages
#BiocManager::install("maEndToEnd")
#BiocManager::install("limma")
#BiocManager::install("oligo")
#BiocManager::install("Biobase")
#BiocManager::install("pd.hugene.1.1.st.v1")
#BiocManager::install('EnhancedVolcano')
options(stringsAsFactors = FALSE)
library(oligo)
library(Biobase)
#library(maEndToEnd)
library(limma)



## adapted from https://www.bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html

raw_data_dir = args[1]
#do I need to duplicate FactorValue and Characteristics?
sdrf_location <- args[2]
fullSDRF <- read.delim(sdrf_location, header = F)
#using "which colnames" beacuse @data has some *redacted* colnames with whitespace
characteristic_col = which(fullSDRF[1,] == args[3])
paired_samples = FALSE
if(length(args) > 3){
	if(substr(as.character(args[4]), 1, 1) == 'T'){
		paired_samples = TRUE
	} else if(substr(as.character(args[4]), 1, 1) == 't'){
		paired_samples = TRUE
	} else if(substr(as.character(args[4]), 1, 1) == '1'){
		paired_samples = TRUE
	}
}
######################################
################# Is data paired?  No. Code below would be for paired data:
################patients_col = match('Characteristics[patient]', fullSDRF[1,]) #match(A, B) returns the index of the first match of 'A value' in 'B vector' or NA for no match.
################paired_samples = TRUE #Just defining the variable
################if(is.na(patients_col) == FALSE){  #there is a column for patients
################	repetitions = table(fullSDRF[-1,patients_col])
################	for(i in repetitions){
################		if(i[[1]] != 2){  #check if each patient appears exactly two times
################			paired_samples = FALSE  #if one fails, its not paired.
################		}
################	}
################} else {
################	paired_samples = FALSE
################}
################print(paste0('paired_samples == ', as.character(paired_samples)))
#####################################

SDRF = read.delim(sdrf_location)[,c(1, 2, 3, 4, characteristic_col)]
colnames(SDRF) = c("file", "ID", "title", "patient", "FactorValue")  #those four first columns are mandatory
rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)
raw_data <- oligo::read.celfiles(filenames = SDRF$file,  verbose = FALSE, phenoData = SDRF)

stopifnot(validObject(raw_data))


Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("ID",
                                     "patient",
                                     "FactorValue")]

exp_raw <- log2(Biobase::exprs(raw_data))


library(ggplot2)

PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                    Tolerance = pData(raw_data)$FactorValue,
                    Individual = pData(raw_data)$patient)

####### Plotting PCA
ggplot(dataGG, aes(PC1, PC2)) +
      geom_point(aes(colour = Tolerance)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  #coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15))
  #+ scale_color_manual(values = c("darkorange2", "dodgerblue4"))



######## Plotting Barplot
###pdf("plot.pdf")   
###oligo::boxplot(raw_data, target = "core", 
###               main = "Boxplot of log2-intensitites for the raw data", las = 2)  
 
# borders might need adjustments


###dev.off()

###########################
# getting the rma normalization

eset_norm  <- oligo::rma(raw_data, target = "core")
expression <- Biobase::exprs(eset_norm)


#### plottng a PCA
###PCA <- prcomp(t(expression), scale = FALSE)
###percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
###sd_ratio <- sqrt(percentVar[2] / percentVar[1])

###dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
###                    FactorValue = 
###                     Biobase::pData(eset_norm)$FactorValue,

###					Individual = Biobase::pData(eset_norm)$patient)

###ggplot(dataGG, aes(PC1, PC2)) +
###      geom_point(aes(colour = FactorValue)) +
###  ggtitle("PCA plot of the calibrated, summarized data") +
###  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
###  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
###  theme(plot.title = element_text(hjust = 0.5)) +
###  coord_fixed(ratio = sd_ratio) #+
###  #scale_shape_manual(values = c(4,15)) + scale_color_manual(values = c("#a00000", "#0000a0", "#909000", "#a0a0a0" , "#a0a0a0" , "#a0a0a0", "#ff0000"))

###ggsave("~/Uploads/myplot2.pdf"
###  , width = 15,
###  height = 15,
###  units = "in")

#########################################
# generate html data quality report

###arrayQualityMetrics(expressionset = raw_data,
###    outdir = tempdir(),
###    force = TRUE, do.logtransform = TRUE,
###    intgroup = c("FactorValue"))


####################################################
# Filtering low intensity probes

medians <- rowMedians(Biobase::exprs(eset_norm))

pdf('histo.pdf')
hist_res <- hist(medians, 100, col = "cornsilk1", freq = FALSE, 
            main = "Histogram of the median intensities", 
            border = "antiquewhite4",
            xlab = "Median intensities")


###dev.off()
#no_of_samples <- 
#  table(paste0(pData(eset_norm)$FactorValue..TOLERANCE., "_", 
#                  pData(eset_norm)$FactorValue..TIME.POINT.))


man_threshold = 4  #type here a manual threshold


samples_cutoff <- 1 #min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(eset_norm), 1,
                           function(x){
                          sum(x > man_threshold) >= samples_cutoff})

table(idx_man_threshold)

exp_manfiltered <- subset(eset_norm, idx_man_threshold)
expression <- Biobase::exprs(exp_manfiltered)

colnames(expression) = pData(raw_data)$patient



########### translating annotation myself

library(dplyr)

organiezed_data = data.frame(expression)
organiezed_data$Probe.Set.ID = rownames(expression)
my_anno_db =  read.csv('MoGene2_081020W_annotation.tsv', sep='\t', stringsAsFactors = F, colClasses = "character")

joined <- left_join(organiezed_data, my_anno_db)


ids_to_exlude <- is.na(joined$Gene.Symbol)
table(ids_to_exlude)
data_frame_final <- joined[(!is.na(joined$Gene.Symbol)) , (names(joined) != 'Probe.Set.ID')]
###Convert to matrix
exp_final = as.matrix(data_frame_final[, (names(data_frame_final) != 'Gene.Symbol')])
rownames(exp_final) = data_frame_final$Gene.Symbol


#########################################################################



####### Another PCA
###PCA <- prcomp(t(expression), scale = FALSE)

###percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
###sd_ratio <- sqrt(percentVar[2] / percentVar[1])

###dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
###                    FactorValue = 
###                     Biobase::pData(exp_manfiltered)$FactorValue,
###					Individual = Biobase::pData(exp_manfiltered)$patient)

###ggplot(dataGG, aes(PC1, PC2)) +
###      geom_point(aes(colour = FactorValue)) +
###  ggtitle("PCA plot of the calibrated, summarized data") +
###  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
###  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
###  theme(plot.title = element_text(hjust = 0.5)) +
###  coord_fixed(ratio = sd_ratio) +
###  scale_shape_manual(values = c(4,15)) +
###  scale_shape_manual(values = c(4,15)) #+ scale_color_manual(values = c("#a00000", "#0000a0", "#909000", "#a0a0a0" , "#a0a0a0" , "#a0a0a0", "#ff0000"))



###ggsave("~/Uploads/myplot2_intens_filtered.pdf"
###  , width = 15,
###  height = 15,
###  units = "in")
  
  
####################3
#probes with multiple mappings
#library(hugene10sttranscriptcluster.db)
##Annotation
#anno <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
#                                  keys = (featureNames(exp_manfiltered)),
#                                  columns = c("SYMBOL", "GENENAME"),
#                                  keytype = "PROBEID")

#anno <- subset(anno, !is.na(SYMBOL))


#library(dplyr)







##########################################333
## Extract data as vectors/data.frame instead of those problematic S4 objects

#########################
#data linear model
#(limma part)



###############################Data not paired
# from limma manual page 37
if(!paired_samples){
	design <- model.matrix(~ 0+factor(SDRF$FactorValue))
	colnames(design) = levels(factor(SDRF$FactorValue))
	#cont.matrix =makeContrasts(Non_Tolerant-Tolerant, levels = design)
	###### The line below was made to avoid typing the content of the SDRF files in the code (unlike commented line above)
	cont.matrix = matrix(c(1, -1), 2, 1, dimnames = list(Levels = levels(factor(SDRF$FactorValue)), 
				Contrasts = paste(levels(factor(SDRF$FactorValue))[[1]], levels(factor(SDRF$FactorValue))[[2]], sep = ' - ')))
	fit <- eBayes(contrasts.fit(lmFit(exp_final, design), cont.matrix))
	table_adj_p <- topTable(fit, adjust="BH", number = Inf)  #Inf to get the whole table

		
}
############################### Paired data
# from limma manual page 44
if(paired_samples){
	designP <- model.matrix(~ 0 + factor(SDRF$FactorValue) + factor(SDRF$patient))
	colnames(designP)[1:2] = c(make.names(levels(factor(SDRF$FactorValue))[[1]]), make.names(levels(factor(SDRF$FactorValue))[[2]]))
	colnames(designP)[3:length(colnames(designP))] = make.names(colnames(designP)[3:length(colnames(designP))])
	cont.matrixP = makeContrasts(paste(colnames(designP)[1], colnames(designP)[2], sep = '-'), levels = designP)
	fitP <- eBayes(contrasts.fit(lmFit(exp_final, designP), cont.matrixP))
	table_adj_p <- topTable(fitP, adjust="BH", number = Inf)  #Inf to get the whole table
}


#plot(table_adj_p[order(table_adj_p$PROBEID),]$AveExpr, table_adj_pP[order(table_adj_pP$PROBEID),]$AveExpr)  
#plot(table_adj_p[order(table_adj_p$PROBEID),]$logFC, table_adj_pP[order(table_adj_pP$PROBEID),]$logFC)
#plot(-1*log10(table_adj_p[order(table_adj_p$PROBEID),]$P.Value), -1*log10(table_adj_pP[order(table_adj_pP$PROBEID),]$P.Value))

final_table_adj_p = table_adj_p[!is.na(table_adj_p$ID),]   # it was SYMBOL on the other datasets


############################ write data
### 
splitted = strsplit(args[2], '.', fixed = T)[[1]]
output_name = paste(paste(splitted[-1*length(splitted)],collapse="."), 'tsv', sep = '.')
write.table(final_table_adj_p, file=output_name, quote=FALSE, sep='\t', row.names = F)


#volcano_names <- ifelse(abs(fit2$coefficients)>=1, 
#                        fit2$genes$SYMBOL, NA)
#volcanoplot(fit2, style = "p-value", highlight = 200, 
#            names = volcano_names,
#            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)



library(EnhancedVolcano)

  EnhancedVolcano(final_table_adj_p,
    lab = final_table_adj_p$ID,  # it was SYMBOL on the other datasets
    x = 'logFC',
    y = 'P.Value',
    pCutoff = 0.01,
    ylim = c(0, 4),
    FCcutoff = 1)
    
