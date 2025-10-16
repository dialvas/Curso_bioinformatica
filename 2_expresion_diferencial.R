##### instalar paquetes de bioconductor
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DESeq2")

#### cargar paquetes 
library(dplyr)
library(DESeq2)
library(gplots)

###### expresion diferencial de genes
#### fuente de los datos https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158854

counts_matrix<-read.table("Data/GSE158854_raw_counts_GRCh38.p13_NCBI.tsv/GSE158854_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
annot<- read.table("Data/anotaciones.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(counts_matrix)

### limpiar y pegar anotaciones
##### vectores para saber a que grupo pertenece la muestra
null <- c("GSM4812255", "GSM4812257", "GSM4812259", "GSM4812260",
          "GSM4812263", "GSM4812265", "GSM4812268")

# Vector con las columnas del grupo Postpartum
postp <- c("GSM4812256", "GSM4812258", "GSM4812261", "GSM4812262",
           "GSM4812264", "GSM4812266", "GSM4812267", "GSM4812269", "GSM4812270")


##### filtrar genes con muy baja expresion- bajos counts

counts_matrix<-counts_matrix%>%
  rowwise() %>%  # para procesar fila por fila
  filter(
    sum(c_across(all_of(null)) >= 10) >= 3 &   # para que al menos 3 columnas en null ≥ 0.5
      sum(c_across(all_of(postp)) >= 10) >= 3   # para que al menos 3 columnas en postp ≥ 0.5
  ) %>%
  ungroup()


exp_matrix_ann<-counts_matrix%>%
  inner_join(annot, by = "GeneID")%>%
  as.data.frame()

##### verificar que los datos de genesymbol son únicos
length(unique(exp_matrix_ann$Symbol))
rm(counts_matrix)
rm(annot)

# Asignar los nombres de las filas con los nombres de los genes
rownames(exp_matrix_ann) <- exp_matrix_ann$Symbol
### elimino columnas que son de anotaciones
exp_matrix_ann$GeneID <- NULL
exp_matrix_ann$Symbol <- NULL
exp_matrix_ann$EnsemblGeneID <- NULL


#### Necesitamos saber QUÉ MUESTRAS pertenecen a QUÉ GRUPO para poder comparar entre ellos.
#### crear dataframe de metadatos
coldata <- data.frame(
  sample = colnames(exp_matrix_ann),
  condition = ifelse(colnames(exp_matrix_ann) %in% null, "nulli", "posp"),
  row.names = colnames(exp_matrix_ann)
)

# Establecer nivel de referencia (posp será la referencia)
coldata$condition <- factor(coldata$condition, levels = c("posp", "nulli"))

#La función DESeqDataSetFromMatrix() crea un objeto de clase DESeqDataSet, 
#que es la estructura de datos fundamental requerida por el paquete DESeq2 para 
#realizar análisis de expresión diferencial en datos de RNA-seq.
# Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = exp_matrix_ann,
  colData = coldata,
  design = ~ condition
)

print(dds)

##Correr en analisis diferencial
dds <- DESeq(dds)
print(sizeFactors(dds))##Extrae y muestra los factores de normalización


# Comparación: nulli vs posp

res <- results(
  dds,
  contrast = c("condition", "nulli", "posp"),
  alpha = 0.05  # FDR threshold  adj p value
)

summary(res)

# Convertir a dataframe
# log2FC > 0 = mayor expresión en nulli
# log2FC < 0 = mayor expresión en posp
res_df <- as.data.frame(res)


# Para visualizaciones (heatmaps, etc.)
normalized_counts <- counts(dds, normalized = TRUE)

sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj)

#Guardar solo genes significativos
sig_norm_counts <- normalized_counts[rownames(sig_genes), ]

# Esto centra cada gen en 0 y muestra desviaciones
heatmap_matrix_scaled <- t(scale(t(sig_norm_counts)))

#### crear heatmap

pdf("results/heatmap.pdf", width = 15, height = 15)
heatmap.2(heatmap_matrix_scaled,col = colorRampPalette(c("blue", "white", "red"))(75), scale="row", 
          key=TRUE, symkey=FALSE,trace="none", cexRow=0.5,margins = c(8, 8) )
dev.off()


#####obtener matriz del heatmap###
### obtener listado de genes para enriquecimiento
### EXPORTAR LISTADO DE GENES DIFERENCIALES -----------------

write.csv(rownames(heatmap_matrix_scaled), file = "results/genes_heatmap.txt", row.names = FALSE, quote = FALSE)

