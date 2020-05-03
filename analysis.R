# Carga de librerías
library(GEOquery)
library(ggplot2)
library(ggrepel)
library(gplots)
library(oligo)
library(Biobase)
library(genefilter)
library(clariomshumantranscriptcluster.db)
library(limma)
library(ReactomePA)
library(reactome.db)
library(org.Hs.eg.db)

# Descarga de los datos desde GEO (ya normalizados)
data_normalized <- getGEO('GSE103157')[[1]]
rownames(pData(data_normalized)) <- pData(data_normalized)$description
colnames(data_normalized) <- rownames(pData(data_normalized))
data_normalized$group <- substr(as.character(data_normalized$description), 
                                1, nchar(as.character(data_normalized$description)) - 3)

# Inspección de los datos normalizados
# boxplot: dado que los los niveles de expresión están normalizados los boxplots de cada muestra son muy similares
boxplot(exprs(data_normalized), las = 2, which = 'all', xaxt = 'n',
        col = c(rep('#efd7ca', 2), rep('#e3b7b1', 2), rep('#945260', 2), rep('#d0ded8', 2), rep('#588b76', 2), rep('#204c39', 2)),
        main = 'Distribución de las intensidades normalizadas de cada muestra')
text(x = 1:12, y = par("usr")[3] - 0.45,
     labels = data_normalized$description, xpd = NA, srt = 20, cex = 0.3)
# PCA: la principal fuente de variabilidad (43'4 %) se corresponde con la condición del nivel de expresión de la proteína CD62L
pca <- prcomp(t(exprs(data_normalized)), scale = FALSE)
loads <- round(pca$sdev^2/sum(pca$sdev^2)*100, 1)
ggplot(data.frame(pca$x), aes(x=PC1, y=PC2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = 'gray70') +
  geom_vline(xintercept = 0, color = 'gray70') +
  geom_point(aes(color = data_normalized$group), alpha = 0.55, size = 3) +
  coord_cartesian(xlim = c(min(pca$x[,1])-5, max(pca$x[,1])+5)) +
  scale_fill_discrete(name = data_normalized$group) +
  geom_text_repel(aes(y = PC2 + 0.25, label = data_normalized$description), segment.size = 0.25, size = 3) +
  labs(x = c(paste('PC1',loads[1],'%')), y=c(paste('PC2',loads[2],'%'))) + 
  ggtitle('Análisis de Componentes Principales') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c('#efd7ca', '#d0ded8', '#e3b7b1', '#588b76', '#945260', '#204c39'))

# Filtraje no específico
# se pasa a trabajar con 18526 genes (de los 24351 iniciales)
annotation(data_normalized) <- 'clariomshumantranscriptcluster.db'
filtered <- nsFilter(data_normalized,
                     require.entrez = TRUE, remove.dupEntrez = TRUE,
                     var.filter = FALSE,
                     feature.exclude = "^AFFX")
filtered$filter.log
data_normalized_filtered <- filtered$eset

# Identificación de los genes diferencialmente expresados
# matriz de diseño
designMatrix <- matrix(c(rep(c(1,0,0,0,0,0), 2),
                         rep(c(0,1,0,0,0,0), 2),
                         rep(c(0,0,1,0,0,0), 2),
                         rep(c(0,0,0,1,0,0), 2),
                         rep(c(0,0,0,0,1,0), 2),
                         rep(c(0,0,0,0,0,1), 2)),
                       ncol = 6, byrow = TRUE)
colnames(designMatrix) <- unique(data_normalized_filtered$group)
rownames(designMatrix) <- data_normalized_filtered$description

# matriz de contrastes
contrastMatrix <- makeContrasts(
  GOODRESPONDERSvsNONRESPONDERS.LOW = GR_CD62L_low_CD4 - NR_CD62L_low_CD4,
  GOODRESPONDERSvsPARCIALRESPONDERS.LOW = GR_CD62L_low_CD4 - IR_CD62L_low_CD4,
  PARCIALRESPONDERSvsNONRESPONDERS.LOW = IR_CD62L_low_CD4 - NR_CD62L_low_CD4,
  HIGHvsLOW = (GR_CD62L_high_CD4 + IR_CD62L_high_CD4 + NR_CD62L_high_CD4) - (GR_CD62L_low_CD4 + IR_CD62L_low_CD4 + NR_CD62L_low_CD4),
  levels = designMatrix
)

# se estima el modelo
fit <- lmFit(data_normalized_filtered, designMatrix)
# se estiman los contrastes
fit.main <- contrasts.fit(fit, contrastMatrix)
# regularización de la varianza utilizando modelos bayesianos empíricos
fit.main <- eBayes(fit.main)

# genes diferencialmente expresados entre las células T CD4+ CD62Llow de los pacientes que responden bien y los que no responden
topTable_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW <- topTable(fit.main, number = nrow(fit.main), coef = 'GOODRESPONDERSvsNONRESPONDERS.LOW', adjust.method = 'fdr')
topTable_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW$SPOT_ID_1 <- NULL
topTableAnnotated_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW <- cbind(topTable_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW[c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')], 
                                                                 select(clariomshumantranscriptcluster.db, 
                                                                        rownames(topTable_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW), c('SYMBOL', 'ENTREZID', 'GENENAME'))[,-1])
head(topTableAnnotated_DEG_GOODRESPONDERSvsNONRESPONDERS.LOW)

# genes diferencialmente expresados entre las células T CD4+ CD62Llow de los pacientes que responden bien y los que lo hacen parcialmente
topTable_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW <- topTable(fit.main, number = nrow(fit.main), coef = 'GOODRESPONDERSvsPARCIALRESPONDERS.LOW', adjust.method = 'fdr')
topTable_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW$SPOT_ID_1 <- NULL
topTableAnnotated_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW <- cbind(topTable_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW[c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')], 
                                                                 select(clariomshumantranscriptcluster.db, 
                                                                        rownames(topTable_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW), c('SYMBOL', 'ENTREZID', 'GENENAME'))[,-1])
head(topTableAnnotated_DEG_GOODRESPONDERSvsPARCIALRESPONDERS.LOW)

# genes diferencialmente expresados entre las células T CD4+ CD62Llow de los pacientes que responden parcialmente y los que no responden
topTable_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW <- topTable(fit.main, number = nrow(fit.main), coef = 'PARCIALRESPONDERSvsNONRESPONDERS.LOW', adjust.method = 'fdr')
topTable_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW$SPOT_ID_1 <- NULL
topTableAnnotated_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW <- cbind(topTable_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW[c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')], 
                                                                     select(clariomshumantranscriptcluster.db, 
                                                                            rownames(topTable_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW), c('SYMBOL', 'ENTREZID', 'GENENAME'))[,-1])
head(topTableAnnotated_DEG_PARCIALRESPONDERSvsNONRESPONDERS.LOW)

# genes diferencialmente expresados entre las células T CD4+ CD62Lhigh y las células T CD4+ CD62Llow
topTable_DEG_HIGHvsLOW <- topTable(fit.main, number = nrow(fit.main), coef = 'HIGHvsLOW', adjust.method = 'fdr')
topTable_DEG_HIGHvsLOW$SPOT_ID_1 <- NULL
topTableAnnotated_DEG_HIGHvsLOW <- cbind(topTable_DEG_HIGHvsLOW[c('logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B')], 
                                         select(clariomshumantranscriptcluster.db, 
                                                rownames(topTable_DEG_HIGHvsLOW), c('SYMBOL', 'ENTREZID', 'GENENAME'))[,-1])
head(topTableAnnotated_DEG_HIGHvsLOW)

# which(topTableAnnotated_DEG_HIGHvsLOW$ENTREZID == 1236) # CCR7
# which(topTableAnnotated_DEG_HIGHvsLOW$ENTREZID == 6932) # TCF7
# which(topTableAnnotated_DEG_HIGHvsLOW$ENTREZID == 940) # CD28

# volcano plot para el contraste 'HIGHvsLOW'
volcanoplot(fit.main, coef = 4, highlight = 20, style = 'p-value',
            names = select(clariomshumantranscriptcluster.db, rownames(fit.main), c('SYMBOL'))$SYMBOL,
            main = paste("Genes expresados diferencialmente", colnames(contrastMatrix)[4], sep="\n"))
abline(v=c(-1,1), h=-log(0.01), col = 'gray70', lty = 2)

# heatmap para el contraste 'HIGHvsLOW'
probesInHeatmap <- rownames(topTableAnnotated_DEG_HIGHvsLOW[topTableAnnotated_DEG_HIGHvsLOW$adj.P.Val < 0.01,])
HMdata <- exprs(data_normalized_filtered)[rownames(exprs(data_normalized_filtered)) %in% probesInHeatmap,]
geneSymbols <- select(clariomshumantranscriptcluster.db, rownames(HMdata), c("SYMBOL"))
SYMBOLS<- geneSymbols$SYMBOL
rownames(HMdata) <- SYMBOLS
heatmap.2(HMdata,
          Colv = FALSE,
          Rowv = TRUE,
          col = bluered(75),
          main = "Genes expresados diferencialmente \n FDR < 0,1, logFC >=1",
          ColSideColors = c(rep('#efd7ca', 2), rep('#e3b7b1', 2), rep('#945260', 2), rep('#d0ded8', 2), rep('#588b76', 2), rep('#204c39', 2)),
          scale = 'row', 
          trace = 'none', 
          density.info = 'none', 
          cexCol = 0.6,
          cexRow = 0.5,
          key = TRUE,
          srtCol = 30)

# Análisis de significación biológica
DEG_HIGHvsLOW <- topTableAnnotated_DEG_HIGHvsLOW[topTableAnnotated_DEG_HIGHvsLOW$adj.P.Val<0.05,]$ENTREZID
enrich.result <- enrichPathway(gene = DEG_HIGHvsLOW,
                               pvalueCutoff = 0.05,
                               readable = T,
                               pAdjustMethod = 'BH',
                               organism = 'human',
                               universe = union(mappedkeys(org.Hs.egGO) , mappedkeys(org.Hs.egPATH)))
enrich.result$Description
