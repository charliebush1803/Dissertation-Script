



library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(ggpubr)
library(tibble)
library(readr)
library(tibble)
 




# Load data
count <- main_derivation_data$raw_counts
info <- main_derivation_data$info

# Inspect
head(count)
head(info)

rownames(count) <- count$ID
countx <- count[, !colnames(count) %in% c("ID", "Symbol")]
countx <- as.matrix(countx)


# Created group variable

group <- factor(info$state)
table(group)


# DGEList

x <- DGEList(counts = countx, group = group)

x <- calcNormFactors(x)

minCpm <- 1
minSample <- 3  

keep <- rowSums(cpm(x) > minCpm) >= minSample
x <- x[keep, , keep.lib.sizes = FALSE]


# MDS Plot of Senescent vs CTRL

col_col <- ifelse(group == "ctr", "darkred", "navy")

plotMDS(x, col = col_col)

# legend figures 

legend("bottomleft",
       legend = c("Control", "Senescent"),
       col = c("darkred", "navy"),
       pch = 16)


table(group)
head(info)
table(info$cell_type)


# MDS Plot of Cell type

plotMDS(x, col = as.numeric(factor(info$batch)))

plotMDS(x, col = as.numeric(factor(info$cell_type)))

cell_type <- factor(info$cell_type)

col_vec <- as.numeric(cell_type)

plotMDS(x, col = col_vec, pch = 16)


legend("bottomleft",
       legend = levels(cell_type),
       col = 1:length(levels(cell_type)),
       pch = 16,
       title = "Cell Type")


#############################################################################

info$cell_type <- factor(info$cell_type)
info$state <- factor(info$state, levels = c("ctr", "senescent"))

design <- model.matrix(~ cell_type + state, data = info)

colnames(design)

v <- voom(x, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = "statesenescent", number = Inf)

head(results)

design <- model.matrix(~ cell_type * state, data = info)
colnames(design)

#############################################################################

# Prepare data

info$cell_type <- factor(info$cell_type)
info$state <- factor(info$state, levels = c("ctr", "senescent"))

table(info$cell_type, info$state)


# create design matrix

design <- model.matrix(~ cell_type + state, data = info)
colnames(design)


# voom + linear model

v <- voom(x, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

# DE results

results <- topTable(fit, coef = "statesenescent", number = Inf, sort.by = "P")

# gene names as a column
results$gene <- rownames(results)

head(results)


#############################################################################


# Volcano plot

results$threshold <- "Not significant"
results$threshold[results$adj.P.Val < 0.05 & results$logFC > 1] <- "Up in senescent"
results$threshold[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down in senescent"

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "Up in senescent" = "red",
    "Down in senescent" = "blue",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Volcano plot: Senescent vs Control",
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value",
    color = ""
  ) +
  theme_minimal(base_size = 14)


#############################################################################



# Clean up DE Results

DE_results2 <- results[!is.na(results$adj.P.Val), ]
top_numbers <- min(50, nrow(DE_results2))
top_heatmap_genes <- rownames(DE_results2)[1:top_numbers]


# Change Ensembl ID to Gene symbols

symbols <- mapIds(
  org.Hs.eg.db,
  keys = top_heatmap_genes,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

symbols[is.na(symbols)] <- top_heatmap_genes[is.na(symbols)]


# Remove Batch effet caused by cell type

mat_corrected <- removeBatchEffect(
  v$E,
  batch = factor(info$cell_type),
  design = model.matrix(~ state, data = info)
)

matrix2 <- mat_corrected[top_heatmap_genes, , drop = FALSE]


# Remove Genes with zero variance
keep_rows <- apply(matrix2, 1, sd, na.rm = TRUE) > 0
matrix2 <- matrix2[keep_rows, , drop = FALSE]

matrix2_scaled <- t(scale(t(matrix2)))

# Row names to Gene symbols
rownames(matrix2_scaled) <- symbols[rownames(matrix2_scaled)]

# OPTIONAL: ensure unique names (important!)
rownames(matrix2_scaled) <- make.unique(rownames(matrix2_scaled))


# annotations of graph

annotation_col <- data.frame(
  State = factor(info$state),
  CellType = factor(info$cell_type)
)

rownames(annotation_col) <- colnames(v$E)
annotation_col <- annotation_col[colnames(matrix2_scaled), , drop = FALSE]


# Heatmap

pheatmap(
  matrix2_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 50 DE genes (senescent vs control, adjusted for cell type)"
)

# Top 50 genes chosen 

top_n <- 50
top_heatmap_genes <- rownames(results)[1:top_n]

matrix2 <- v$E[top_heatmap_genes, ]

# Z-score by gene (row scaling)
matrix2_scaled <- t(scale(t(matrix2)))

# Annotation for samples
annotation_col <- data.frame(
  State = info$state,
  CellType = info$cell_type
)
rownames(annotation_col) <- rownames(info)

# Heatmap
pheatmap(
  matrix2_scaled,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = paste("Top", top_n, "DE genes: Senescent vs Control")
)






sig_genes <- results[results$adj.P.Val < 0.05, ]
genes <- rownames(sig_genes)

gene_entrez <- mapIds(
  org.Hs.eg.db,
  keys = genes,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

gene_entrez <- na.omit(gene_entrez)

ego <- enrichGO(
  gene = gene_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",              # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


dotplot(ego, showCategory = 15) +
  ggtitle("GO enrichment: senescence-associated genes")

barplot(ego, showCategory = 15)


## cellular senescence signature 

sig_strong <- results[
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
]

up_genes <- rownames(sig_strong[sig_strong$logFC > 0, ])
down_genes <- rownames(sig_strong[sig_strong$logFC < 0, ])

symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(sig_strong),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)

sig_strong$symbol <- symbols

signature_genes <- unique(na.omit(sig_strong$symbol))
length(signature_genes)

validation_signature <- c("CDKN1A", "CDKN2A" , "IL6", "IL8", "CXCL10", "MMP3", "SERPINE1", "MKI67", "PCNA", "CCNBI", "CDK1" )

overlap <- intersect(signature_genes, validation_signature)

length(overlap)
overlap

fisher_matrix <- matrix(c(
  length(overlap),
  length(signature_genes) - length(overlap),
  length(validation_signature) - length(overlap),
  20000  # approx total genes
), nrow = 2)

fisher.test(fisher_matrix)


################################################################


# Validation dataset 


data <- validation_data1_GSE164012

str(data)

countV <- data$raw_counts
infoV  <- data$info
table(infoV$group)
group <- factor(infoV$group)
table(group)

group <- factor(infoV[colnames(countV), "group"])
levels(group) <- gsub("/", "_", levels(group))


table(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)


head(colnames(countV))
head(rownames(infoV))
head(infoV$GSM)
head(infoV$title)

sum(colnames(countV) %in% rownames(infoV))

idx <- match(colnames(countV), rownames(infoV))
group <- factor(infoV$group[idx])
levels(group) <- gsub("/", "_", levels(group))
table(group)

colnames(countV)
group

#count matrix to data
idx <- match(colnames(countV), rownames(infoV))

idx
colnames(countV)[is.na(idx)]

symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(countV),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# DGE

dge <- DGEList(counts = countV2)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)



contrast_bleo_ctrl <- makeContrasts(BLEO_vs_CTRL = BLEO - CTRL, levels = design)
contrast_bleo_pcc1 <- makeContrasts(BLEO_vs_PCC1 = BLEO - BLEO_PCC1, levels = design)

lrt_bleo_ctrl <- glmLRT(fit, contrast = contrast_bleo_ctrl)
lrt_bleo_pcc1 <- glmLRT(fit, contrast = contrast_bleo_pcc1)

# Extract DE results
de_bleo_ctrl <- topTags(lrt_bleo_ctrl, n = Inf)$table
de_bleo_pcc1 <- topTags(lrt_bleo_pcc1, n = Inf)$table

# Add gene column (using rownames from your matrix)
de_bleo_ctrl$gene <- rownames(de_bleo_ctrl)
de_bleo_pcc1$gene <- rownames(de_bleo_pcc1)


de_bleo_ctrl <- topTags(lrt_bleo_ctrl, n = Inf)$table %>%
  rownames_to_column("gene")

de_bleo_pcc1 <- topTags(lrt_bleo_pcc1, n = Inf)$table %>%
  rownames_to_column("gene")


sig_up <- intersect(
  de_bleo_ctrl %>% filter(FDR < 0.05, logFC > 1) %>% pull(gene),
  de_bleo_pcc1 %>% filter(FDR < 0.05, logFC > 1) %>% pull(gene)
)

sig_down <- intersect(
  de_bleo_ctrl %>% filter(FDR < 0.05, logFC < -1) %>% pull(gene),
  de_bleo_pcc1 %>% filter(FDR < 0.05, logFC < -1) %>% pull(gene)
)

length(sig_up)
length(sig_down)


# Senescence Score 


senescence_score <- function(expr, up_genes, down_genes = NULL) {
  
  up_present <- intersect(up_genes, rownames(expr))
  up_score <- colMeans(expr[up_present, , drop = FALSE], na.rm = TRUE)
  
  if (!is.null(down_genes)) {
    down_present <- intersect(down_genes, rownames(expr))
    down_score <- colMeans(expr[down_present, , drop = FALSE], na.rm = TRUE)
    score <- up_score - down_score
  } else {
    score <- up_score
  }
  
  return(score)

  
  dge <- DGEList(counts = countV2)
  dge <- calcNormFactors(dge)
  
  expr_norm <- cpm(dge, log = TRUE, prior.count = 1)
  
  group <- factor(group, levels = c("CTRL", "BLEO_PCC1", "BLEO"))
  
  up_present <- intersect(sig_up, rownames(expr_norm))
  down_present <- intersect(sig_down, rownames(expr_norm))
  
  cat("Up genes present:", length(up_present), "\n")
  cat("Down genes present:", length(down_present), "\n")
  
  if (length(up_present) < 3) stop("Too few up genes present to build a stable score.")
  if (length(down_present) < 3) warning("Very few down genes present; score may be less stable.")
  

  z_expr <- t(scale(t(expr_norm)))
  

  z_expr <- z_expr[complete.cases(z_expr), , drop = FALSE]
  
  up_present <- intersect(up_present, rownames(z_expr))
  down_present <- intersect(down_present, rownames(z_expr))
  

  up_score <- colMeans(z_expr[up_present, , drop = FALSE], na.rm = TRUE)
  
  if (length(down_present) > 0) {
    down_score <- colMeans(z_expr[down_present, , drop = FALSE], na.rm = TRUE)
    senescence_score <- up_score - down_score
  } else {
    senescence_score <- up_score
  }
  
 
  group_means <- tapply(senescence_score, group, mean, na.rm = TRUE)
  
  if (group_means["BLEO"] < group_means["CTRL"]) {
    senescence_score <- -senescence_score
    group_means <- tapply(senescence_score, group, mean, na.rm = TRUE)
  }
  
  print(group_means)
  
  ## 7. Save the final score table
  senescence_scores <- data.frame(
    Sample = colnames(expr_norm),
    Group = group,
    SenescenceScore = as.numeric(senescence_score)
  )
  
  aggregate(SenescenceScore ~ Group, data = senescence_scores, mean)
  
  
  senescence_scores$Group <- factor(
    senescence_scores$Group,
    levels = c("CTRL", "BLEO_PCC1", "BLEO")
  )
  
  comparisons <- list(
    c("CTRL", "BLEO"),
    c("CTRL", "BLEO_PCC1"),
    c("BLEO_PCC1", "BLEO")
  )
  
  p <- ggplot(senescence_scores, aes(x = Group, y = SenescenceScore, fill = Group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.08, size = 2.5, alpha = 0.9) +
    stat_compare_means(
      comparisons = comparisons,
      method = "t.test",
      label = "p.format",
      hide.ns = FALSE,
      step.increase = 0.1
    ) +
    labs(
      title = "Senescence score by group",
      x = NULL,
      y = "Senescence score"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(face = "bold"),
      legend.position = "none"
    )
  
  print(p)

  
  ctrl_mean <- mean(senescence_scores$SenescenceScore[senescence_scores$Group == "CTRL"])
  
  senescence_scores$Score_centered <- senescence_scores$SenescenceScore - ctrl_mean
  
  aes(y = Score_centered)
  
  ylab("Senescence score (z-score)")
  
  comparisons <- list(
    c("CTRL", "BLEO"),
    c("BLEO_PCC1", "BLEO")


# Signature matrix in a Table

    signature_matrix <- data.frame(
      Gene = rownames(expr_norm),
      Senescent = rowMeans(expr_norm[, group == "BLEO"]),
      Non_senescent = rowMeans(expr_norm[, group != "BLEO"])
    )
    
    write.table(signature_matrix,
                "senescence_signature_matrix.txt",
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    
    
    head(signature_matrix)


##########################################################
    
#  Upload to CIBERSORT
    
    
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    
packages_cran <- c("dplyr", "readr", "tibble")
packages_bioc <- c("edgeR", "org.Hs.eg.db", "AnnotationDbi")
    
for (pkg in packages_cran) {
if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    }
for (pkg in packages_bioc) {
if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    
    
    
    
    

data3 <- main_derivation_data
count3 <- data$raw_counts
info3  <- data$info

# Optional: use existing CPM if present
normcpm <- data$normalised_cpm

## ---------------------------
## 3. Basic checks
## ---------------------------
cat("Count matrix dimensions:", dim(count3), "\n")
cat("Info dimensions:", dim(info3), "\n")
cat("First few sample names in count matrix:\n")
print(head(colnames(count3)))
cat("First few rows of metadata:\n")
print(head(info3))

# Check sample grouping
if (!"state" %in% colnames(info3)) {
    
    
  table(info$state)
  
  if (!all(colnames(count) %in% rownames(info))) {

    if (all(colnames(count) %in% rownames(info))) {
      info <- info[colnames(count), , drop = FALSE]
    }
    
    rownames(count) <- sub("\\..*$", "", rownames(count))   
    
    gene_symbols <- mapIds(
      org.Hs.eg.db,
      keys = rownames(count),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )    
    
  
    
    rownames(count) <- count$ID

count$ID <- NULL

rownames(count) <- sub("\\..*$", "", rownames(count))

head(rownames(count))

library(org.Hs.eg.db)
library(AnnotationDbi)

symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(count),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

keep <- !is.na(symbols) & symbols != ""
count <- count[keep, , drop = FALSE]
rownames(count) <- symbols[keep]

count <- count[!duplicated(rownames(count)), , drop = FALSE]

head(rownames(count))
######################################################################################################################################################

annot_df <- data.frame(
  ENSEMBL = rownames(count),
  SYMBOL = unname(gene_symbols),
  stringsAsFactors = FALSE
)

annot_df <- annot_df %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

count_mapped <- count[annot_df$ENSEMBL, , drop = FALSE]
rownames(count_mapped) <- annot_df$SYMBOL  
    
    
    
    
    
########################################################################################    

write.table(sig_up, "senescence_signature_up.txt",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
write.table(sig_down, "senescence_signature_down.txt",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    signature_matrix <- data.frame(
      GeneSymbol = rownames(expr_norm),
      Senescent = rowMeans(expr_norm[, group == "BLEO", drop = FALSE]),
      Non_senescent = rowMeans(expr_norm[, group %in% c("CTRL", "BLEO_PCC1"), drop = FALSE]),
      check.names = FALSE
    )
    
    sig_all <- unique(c(sig_up, sig_down))
    signature_matrix <- signature_matrix[signature_matrix$GeneSymbol %in% sig_all, ]
    
    signature_matrix <- signature_matrix[!duplicated(signature_matrix$GeneSymbol), ]
    
    write.table(
      signature_matrix,
      file = "senescence_signature_matrix.txt",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )


    
    expr <- as.data.frame(expr)
    
   
    expr$GeneSymbol <- rownames(expr)
    
    expr <- expr[, c("GeneSymbol", setdiff(colnames(expr), "GeneSymbol"))]
    
    expr <- expr[!is.na(expr$GeneSymbol) & expr$GeneSymbol != "", ]
    
    library(dplyr)
    
    expr_clean <- expr %>%
      group_by(GeneSymbol) %>%
      summarise(across(everything(), median, na.rm = TRUE), .groups = "drop")
    
    write.table(
      expr_clean,
      file = "CIBERSORTx_mixture.txt",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    
    # remaining columns = samples
    expr_out <- data.frame(GeneSymbol = rownames(expr), expr, check.names = FALSE)
    
    write.table(
      expr_out,
      file = "CIBERSORTx_mixture_senescence_ready.txt",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.csv(
      senescence_score_df,
      file = "Senescence_sample_scores.csv",
      row.names = FALSE
    )
    
    saveRDS(expr, file = "CIBERSORTx_mixture_senescence_ready.rds")
    
    cat("Files written:\n")
    cat("- CIBERSORTx_mixture_senescence_ready.txt\n")
    cat("- Senescence_sample_scores.csv\n")
    cat("- CIBERSORTx_mixture_senescence_ready.rds\n")
    
    ## ---------------------------
    ## 10. Optional: build a phenotype table
    ## ---------------------------
    if ("state" %in% colnames(info) && all(colnames(expr) %in% rownames(info))) {
      phenotype_table <- data.frame(
        Sample = colnames(expr),
        State = info[colnames(expr), "state"],
        SenescenceScore = senescence_score[colnames(expr)],
        row.names = NULL
      )
      
      write.csv(phenotype_table, "Sample_phenotype_table.csv", row.names = FALSE)
      cat("- Sample_phenotype_table.csv\n")
    }



