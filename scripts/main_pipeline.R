# Set up ----

# Install CRAN packages
install.packages(c("colorspace", "dplyr", "ggraph", "ggplot2", "igraph", "pheatmap", "readr", "scales", "scatterpie", "shadowtext", "stringr", "tidyr", "viridis"))

# Install BiocManager packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.22")

BiocManager::install(c("clusterProfiler", "ComplexHeatmap", "DEqMS", "EnhancedVolcano", "enrichplot", "limma", "msImpute", "org.Hs.eg.db", "pathview", "PRONE", "RCy3", "rrvgo", "STRINGdb"))

# Data handeling
handeling <- c("readr", "tidyr", "dplyr", "stringr")
lapply(handeling, library, character.only = TRUE)

# Visualization
visualization <- c("ggplot2", "grid", "pheatmap", "ComplexHeatmap", "EnhancedVolcano", "scales", "shadowtext", "colorspace", "RColorBrewer", "viridis", "grid")
lapply(visualization, library, character.only = TRUE)

# Quality Control and Preprocessing
QC_preprocessing <- c("PRONE", "msImpute")
lapply(QC_preprocessing, library, character.only = TRUE)

# Differential expression analysis
dea<- c("limma", "DEqMS")
lapply(dea, library, character.only = TRUE)

# Enrichment analysis
ea<- c("clusterProfiler", "org.Hs.eg.db", "rrvgo","enrichplot")
lapply(ea, library, character.only = TRUE)

# Pathway visualization
library(pathview)

# Network analysis
na<- c("STRINGdb", "igraph", "ggraph", "scatterpie", "RCy3")
lapply(na, library, character.only = TRUE)


# Data loading ----

MQ_data <- read_tsv("proteinGroups.txt")
dim(MQ_data)
str(MQ_data)


md <- read_tsv("metadata.txt")
dim(md)
str(md)


# Quality Control ----

data <- MQ_data %>%
  dplyr::select(
    # Identifiers
    `Majority protein IDs`,
    `Gene names`,
    
    # Quantification
    matches("^Reporter intensity corrected .*Pat_tbc_pool_"),
    
    # QC / filtering
    `Peptide counts (razor+unique)`,
    `Only identified by site`,
    `Reverse`,
    `Potential contaminant`
  )

data <- data %>%
  dplyr::mutate(
    # Representative protein ID
    protein_id = str_split_fixed(`Majority protein IDs`, ";", 2)[, 1],
    
    # Representative gene symbol
    gene = str_split_fixed(`Gene names`, ";", 2)[, 1],
    
    # Representative peptide count
    n_peptides = as.integer(
      str_split_fixed(`Peptide counts (razor+unique)`, ";", 2)[, 1]
    )
  ) %>%
  # Drop the original columns to avoid ambiguity
  dplyr::select(-`Majority protein IDs`, -`Gene names`, -`Peptide counts (razor+unique)`)

anyDuplicated(data$protein_id) # Should output 0 if there are no duplicates


data_qc <- data %>%
  
  # Remove contaminants and reverse hits
  dplyr::filter(
    is.na(`Only identified by site`) | `Only identified by site` != "+",
    is.na(`Reverse`) | `Reverse` != "+",
    is.na(`Potential contaminant`) | `Potential contaminant` != "+"
  ) %>%
  
  # (Optional) Remove proteins identified only by 1 peptide
  dplyr::filter(n_peptides >= 2) %>%
  
  # Drop used columns
  dplyr::select(-`Only identified by site`, -`Reverse`, -`Potential contaminant`)


data_qc <- data_qc %>% mutate(across(starts_with("Reporter intensity corrected"), ~ na_if(., 0)))


# Select quantification columns
sample_cols <- grep("^Reporter intensity corrected", names(data_qc))

# Drop features with to many NAs
na_counts <- rowSums(is.na(data_qc[, sample_cols]))
proteins_to_keep <- na_counts < (0.3 * length(sample_cols)) 
data_qc <- data_qc[proteins_to_keep, ]

# Drop samples with to many NAs
data_qc <- data_qc %>%
  dplyr::select(
    -all_of(sample_cols),
    all_of(sample_cols[colMeans(is.na(dplyr::select(., all_of(sample_cols)))) < 0.3])
  )

# Filter the metadata to keep same samples
md_qc <- md %>%
  dplyr::filter(Column %in% intersect(md$Column, colnames(data_qc))) %>%
  dplyr::mutate(Column = factor(Column, levels = colnames(data_qc))) %>%
  dplyr::arrange(Column) %>%
  dplyr::mutate(Column = as.character(Column))


# Normalization ----


se <- load_data(data_qc, # Expression matrix
                md_qc, # Metadata
                protein_column = "protein_id", 
                gene_column = "gene", 
                ref_samples = md[md$Condition == "ref",]$Column, 
                batch_column = "Batch", 
                condition_column = "Condition", 
                label_column = "Label")

se_norm <- normalize_se(se, c("Median"  # Simple sample shifting
                              ,"Quantile"# Sample-to-reference
                              ,"LoessF"  # Sample-to-reference
                              ,"VSN"     # Variance stabilizing
                              ,"RobNorm" # Model-based
))

plot_boxplots(se_norm, 
              ain = NULL, 
              color_by = NULL, 
              label_by = NULL, 
              ncol = 3, 
              facet_norm = TRUE) + 
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")

plot_intragroup_PMAD(se_norm, 
                     ain = NULL, 
                     condition = NULL, 
                     diff = TRUE) +
  ggplot2::ylab("Intragroup variation compared to log2") +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")


# Select the desired assay
abund <- as.matrix(assay(se_norm, "log2"))
abund <- abund[complete.cases(abund), ]

# Calculate PCs
pca_res <- prcomp(t(abund), center = TRUE, scale. = TRUE)

# Extract PCA results
scores <- as.data.frame(pca_res$x)
scores$Column <- rownames(scores)
scores <- merge(scores, md, by = "Column", all.x = TRUE)

# Variance explained for axis labels
ve <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

# Plotting code
ggplot(scores, aes(x = PC1, y = PC2, color = Batch, shape = Condition)) +
  geom_point(size = 5) +
  labs(
    title = " log2",
    x = paste0("PC1 (", round(100 * ve[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * ve[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 10),
    legend.position = "right"
  )

se_norm <- normalize_se(se, c("RobNorm", "IRS_on_RobNorm", "limBE_on_RobNorm"), combination_pattern = "_on_")

# Select the desired assay
abund <- as.matrix(assay(se_norm, "limBE_on_RobNorm"))
abund <- abund[complete.cases(abund), ]

# Calculate PCs
pca_res <- prcomp(t(abund), center = TRUE, scale. = TRUE)

# Extract PCA results
scores <- as.data.frame(pca_res$x)
scores$Column <- rownames(scores)
scores <- merge(scores, md, by = "Column", all.x = TRUE)

# Variance explained for axis labels
ve <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

# Plotting code
ggplot(scores, aes(x = PC1, y = PC2, color = Batch, shape = Condition)) +
  geom_point(size = 5) +
  labs(
    title = "limBE on RobNorm",
    x = paste0("PC1 (", round(100 * ve[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * ve[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text  = element_text(size = 10),
    legend.position = "right"
  )


data_norm <- assay(se_norm, "limBE_on_RobNorm")

# Rename the columns with sample label
label_map <- setNames(md$Label, md$Column)
new_labels <- label_map[colnames(data_norm)]
colnames(data_norm) <- new_labels

data_norm <- as.matrix(data_norm)

# Recover feature ids
rownames(data_norm) <- rowData(se_norm)$Protein.IDs


# Calculate correlations
cor_mat <- cor(data_norm, use = "pairwise.complete.obs", method = "pearson")

# Plot correlation matrix
pheatmap::pheatmap(
  cor_mat,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation"
)

# NA imputation ----

md_qc<-as.data.frame(md_qc)
batch <- factor(md_qc[, "Batch"])
condition <- factor(md_qc[, "Condition"])

hdp <- selectFeatures(data_norm, method = "ebm", group = condition)

# Construct matrix M to capture missing entries
M <- ifelse(is.na(data_norm),1,0)
M <- M[hdp$msImpute_feature,]

# Build the annotation
ann_col <- md %>%
  filter(Label %in% colnames(M)) %>%
  dplyr::select(Label, Batch, Condition) %>%
  distinct() %>%
  as.data.frame()

rownames(ann_col) <- ann_col$Label
ann_col$Label <- NULL

# Ensure factors
ann_col$Batch <- factor(ann_col$Batch)
ann_col$Condition <- factor(ann_col$Condition)

# Build annotation colors
cond_lvls <- levels(ann_col$Condition)
batch_lvls <- levels(ann_col$Batch)

ann_colors <- list(
  Condition = setNames(brewer.pal(max(3, length(cond_lvls)), "Set1")[seq_along(cond_lvls)], cond_lvls),
  Batch     = setNames(brewer.pal(max(3, length(batch_lvls)), "Set2")[seq_along(batch_lvls)], batch_lvls)
)

# Plotting code
ha <- HeatmapAnnotation(
  df  = ann_col,
  col = ann_colors,
  annotation_name_side = "right",
  gp = gpar(col = "white", lwd = 0.5)
)

hm <- Heatmap(M,
              column_title = "NA pattern, columns ordered by dropout similarity",
              name = "Pressence",
              col = c("#FFEFDB", "#8FBC8F"),
              show_row_names = FALSE,
              show_column_names = FALSE,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              top_annotation = ha,
              row_names_gp =  gpar(fontsize = 7),
              column_names_gp = gpar(fontsize = 8),
              heatmap_legend_param = list(#direction = "horizontal",
                labels = c("missing", "observed"),
                legend_width = unit(6, "cm"))#,
              #rect_gp = gpar(col = NA, lwd = 0.5)
)

hmd <- draw(hm,
            heatmap_legend_side = "right",
            annotation_legend_side = "right",
            merge_legends = TRUE
)


msImpute_mar <- msImpute(data_norm, method = "v2") # method = "v2" is for MAR assumptions


# Sort the plots in the display
par(mfrow=c(1,2))

# Generate the CV plots
pcv <- plotCV2(data_norm, main = "data")
pcv <- plotCV2(msImpute_mar, main = "msImpute v2")


data <- msImpute_mar ## IMPORTANT


# Diferential Expression Analysis ----

md_qc$Condition <- factor(md_qc$Condition)

design <- model.matrix(~ 0 + Condition, data = md_qc) # Batch can be modeled if not already taken care in the normalization step

colnames(design) <- sub("^Condition", "", colnames(design)) # nicer names (optional)

fit <- lmFit(data, design)

# Get condition combinations
conds <- colnames(design)
pairs <- combn(conds, 2)

# Create contrasts
contrast_names <- apply(pairs, 2, \(x) paste0(x[1], "_vs_", x[2])) # B_vs_A
contrast_formulas <- apply(pairs, 2, \(x) paste0(x[1], " - ", x[2])) # B - A

Cmat <- makeContrasts(contrasts = contrast_formulas, levels = design)
colnames(Cmat) <- contrast_names

fit2 <- contrasts.fit(fit, Cmat)
fit2 <- eBayes(fit2)

res <- topTable(fit2
                , coef = "HC_vs_PTB" # Comparison of interest
                , adjust.method="BH"
                , sort.by="P"
                , number=Inf)

sig_counts <- data.frame(
  Contrast = colnames(fit2),
  N_sig = sapply(colnames(fit2), function(cn) {
    tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 0, na.rm = TRUE)
  })
)

# Get peptide column
pep_count <- setNames(data_qc$n_peptides, data_qc$protein_id)

# Add the peptide count column
fit2$count <- pep_count[rownames(fit2$coefficients)]

fit3 <- spectraCounteBayes(fit2)

res_deqms <- outputResult(fit3, coef_col = "HC_vs_PTB")

res_deqms$gene <- rowData(se_norm)$Gene.Names[
  match(rownames(res_deqms), rowData(se_norm)$Protein.IDs)
] 

sig_counts_deqms <- data.frame(
  Contrast = colnames(fit2$coefficients),
  N_sig = sapply(colnames(fit2$coefficients), function(cn) {
    tt <- outputResult(fit3, coef_col = cn)   # DEqMS output table for that contrast
    sum(tt$sca.adj.pval < 0.05 & abs(tt$logFC) > 0, na.rm = TRUE)
  })
)


EnhancedVolcano(res_deqms, # DEA results
                x = "logFC", # FC column
                y = "sca.adj.pval", # Adjusted P column
                lab = res_deqms$gene, # Gene symbols for visualization
                labSize = 5, 
                pCutoff = 5e-2, # P cut-off
                FCcutoff = 1, # FC cut-off
                title = "HC vs. PTB",
                subtitle = bquote(italic("DEqMS"))) +
  scale_x_continuous(limits = c(-8,  8), breaks = seq(-8, 8, by = 2)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2))


# Build matrix
top_proteins_deqms <- rownames(res_deqms[abs(res_deqms$logFC) > 1, ])[ order(res_deqms$sca.adj.pval[abs(res_deqms$logFC) > 1]) ][1:20]
cols_keep <- md_qc$Label[md_qc$Condition %in% c("HC", "PTB")]
mat_top_deqms <- data[top_proteins_deqms, cols_keep, drop = FALSE]
rownames(mat_top_deqms) <- res_deqms$gene[match(rownames(mat_top_deqms), rownames(res_deqms))]

# Build annotation
ann_col <- md %>%
  filter(Label %in% colnames(mat_top_deqms)) %>%
  dplyr::select(Label, Batch, Condition) %>%
  as.data.frame()
rownames(ann_col) <- ann_col$Label
ann_col$Label <- NULL
ann_col$Batch <- factor(ann_col$Batch)
ann_col$Condition <- factor(ann_col$Condition)

# Build annotation colors
cond_lvls <- levels(ann_col$Condition)
batch_lvls <- levels(ann_col$Batch)
ann_colors <- list(
  Condition = setNames(brewer.pal(max(3, length(cond_lvls)), "Set1")[seq_along(cond_lvls)], cond_lvls),
  Batch     = setNames(brewer.pal(max(3, length(batch_lvls)), "Set2")[seq_along(batch_lvls)], batch_lvls)
)

# Plotting code
top_ha <- HeatmapAnnotation(
  df  = ann_col,
  col = ann_colors,
  annotation_name_side = "right",
  gp = gpar(col = "white", lwd = 0.5)
)

ht <- Heatmap(
  mat_top_deqms,
  name = "Intensity",
  top_annotation = top_ha,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 11),
  col = viridis(100),
  rect_gp = gpar(col = "white", lwd = 0.5)
)

htt<-draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE
)

# Over Representation Analysis ----

res_sig <- res_deqms %>% dplyr::filter(sca.adj.pval <= 0.05, abs(logFC) >= 1)
all_proteins <- unique(data_qc$gene) 
sig_proteins <- unique(res_sig$gene)

go_enrich <- enrichGO(gene = sig_proteins,
                      universe = all_proteins,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.35,
                      qvalueCutoff = 0.33
)

# Inspect the first couple of terms
head(go_enrich@result[,c("ID","Description","GeneRatio","BgRatio","p.adjust")])

# Create dataframe
enrich_data <- as.data.frame(go_enrich)

dotplot(go_enrich
        , x = "RichFactor"
        , showCategory = 10)

cnetplot(go_enrich
         , showCategory = 2
         , node_label = "none"
         , size_item=0.6
         , size_category=1.2
         , color_category="#97c497"
         , color_item="#c4c4c4") +
  geom_cnet_label(node_label = "category", size = 5) +
  geom_cnet_label(node_label = "item", size = 4)

orgdb="org.Hs.eg.db"
ont="BP"
keytype = "SYMBOL"

# Calculate similarity matrix
simMatrix <- calculateSimMatrix(
  enrich_data$ID,
  orgdb,
  keytype = "SYMBOL",
  ont = "BP",
  method = "Rel",
  semdata = GOSemSim::godata(annoDb = orgdb, ont = ont, keytype = keytype))

# Calculate ranking metric
enrich_data$logp <- (-log10(enrich_data$p.adjust))
enrichment_scores <- setNames(enrich_data$logp, enrich_data$ID)

# Select and group under the parent term
reducedTerms <- reduceSimMatrix(simMatrix,
                                enrichment_scores,
                                threshold=0.8, # Similarity threshold
                                orgdb="org.Hs.eg.db")

reduced_heatmap <- heatmapPlot(simMatrix,
                               reducedTerms,
                               annotateParent=TRUE,
                               annotationLabel="parentTerm",
                               fontsize=6)

# GSEA ----


# Get the desired columns from the DEA result
gene_abundance <- res_deqms %>%
  dplyr::select(logFC) %>%
  dplyr::arrange(desc(logFC)) %>%
  tidyr::drop_na()

# Create the named vector
geneList <- setNames(
  gene_abundance$logFC,
  rownames(gene_abundance)
)

gsea_kegg <- gseKEGG(geneList = geneList, 
                     organism = "hsa", 
                     keyType = "uniprot", 
                     pAdjustMethod = "BH", 
                     pvalueCutoff = 1, # Cutoff set to one to keep all results
                     verbose = FALSE)

gsea_kegg2 <- arrange(gsea_kegg, desc(abs(NES)))
# Plot top 3 gene sets
gseaplot2(gsea_kegg2, 1:3, pvalue_table=T, base_size=14)

# Pathway Vizualization ----

gene_abundance <- res_deqms %>%
  dplyr::select(gene, logFC) %>%
  dplyr::arrange(desc(logFC)) %>%
  tidyr::drop_na()

geneList <- setNames(
  gene_abundance$logFC,
  gene_abundance$gene
)

pv_out <- pathview(gene.data  = geneList,
                   gene.idtype = "SYMBOL",
                   pathway.id = "hsa04518", # pathway of interest from the GSEA analysis
                   species    = "hsa",
                   out.suffix = "test1",
                   limit = list(gene = 6),
                   both.dirs = list(gene = TRUE),
                   low = list(gene = "#2166ac"),
                   mid = list(gene = "#f7f7f7"),
                   high = list(gene = "#b2182b"))

# Network Analysis ----

string_db <- STRINGdb$new(version="11.5"
                          , species=9606
                          , score_threshold=700) # Threshod to filter only top confidence interactions and reduce network cluttering

map_df <- data.frame(UNIPROT=rownames(res_sig))
mapped  <- string_db$map(map_df
                         , "UNIPROT"
                         , removeUnmappedRows = TRUE)
hits <- unique(mapped$STRING_id)

edges <- string_db$get_interactions(hits) # Retrieve the interaction network for these hits

edges$weight <- edges$combined_score / 1000


g <- graph_from_data_frame(
  d = unique(edges[, c("from","to","weight")]),
  directed = FALSE)


# Get variables to add as atributes
node_attr <- res_sig %>%
  tibble::rownames_to_column("UNIPROT") %>%
  transmute(
    UNIPROT,
    gene = gene,
    logFC = logFC,
    fdr = adj.P.Val,
    neglog10fdr = -log10(adj.P.Val)
  )

node_attr_string <- mapped %>%
  left_join(node_attr, by = "UNIPROT")

# Add attributes
m <- match(V(g)$name, node_attr_string$STRING_id)

V(g)$UNIPROT <- node_attr_string$UNIPROT[m]
V(g)$gene <- node_attr_string$gene[m]
V(g)$logFC <- node_attr_string$logFC[m]
V(g)$fdr <- node_attr_string$fdr[m]
V(g)$neglog10fdr <- node_attr_string$neglog10fdr[m]

V(g)$label <- V(g)$gene


# Topology metrics
V(g)$degree <- igraph::degree(g, mode="all")
V(g)$betweenness <- betweenness(g, directed=FALSE, normalized=TRUE)
V(g)$closeness <- closeness(g)
V(g)$eigen <- eigen_centrality(g)$vector

# Proteins with the highest degree are identified as potential hub proteins
top_hubs <- V(g)$name[order(V(g)$degree, decreasing = TRUE)][seq_len(min(20, igraph::vcount(g)))]


g <- igraph::simplify(g, remove.multiple=TRUE, remove.loops=TRUE)
com <- cluster_louvain(g, weights = E(g)$weight)
V(g)$module <- membership(com)


set.seed(27)

nw <- 
  ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = weight), colour = "gray", show.legend = FALSE) +
  geom_node_point(
    aes(fill = logFC, size = degree),
    shape = 21, colour = "black", stroke = 0.6, alpha = 0.95) +
  scale_size_continuous(range = c(4, 15), name = "Degree") +
  scale_fill_gradient2(
    low = "blue", mid = "grey90", high = "red",
    midpoint = 0, name = "log2(FC)"
  ) +
  guides(
    color = "none",
    fill = guide_colorbar(order = 1),
    size  = guide_legend(order = 2)
  ) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right",
        legend.box = "vertical") +
  ggtitle("STRING PPI colored by log2(FC)") +
  geom_shadowtext(
    #data = function(d) d[d$name %in% top_hubs, ],
    aes(x = x, y = y, label = label),
    colour   = "black",
    bg.colour = "white",
    bg.r      = 0.1,
    size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) + coord_fixed()
nw


set.seed(27)

mods <- V(g)$module
mod_ids <- sort(unique(mods))
module_pal <- setNames(hue_pal(l = 65, c = 100)(length(mod_ids)), mod_ids)
module_pal_fill  <- lighten(module_pal, amount = 0.30)
module_pal_border <- darken( module_pal, amount = 0.30)

nw_mod <- 
  ggraph(g, layout = "fr") +
  geom_edge_link(aes(alpha = weight), colour = "grey75", show.legend = FALSE) +
  geom_node_point(aes(size = degree, fill = factor(module), color = factor(module)),
                  shape = 21, stroke = 0.6, alpha = 0.95) +
  scale_size(range = c(4, 15), name = "Degree") +
  scale_fill_manual(values = module_pal_fill,   name = "Module") +
  scale_color_manual(values = module_pal_border, guide = "none") +
  guides(
    color = "none",
    fill = guide_legend(order = 1, override.aes = list(size = 6)),
    size  = guide_legend(order = 2)
  ) +
  theme_void() +
  theme(plot.margin = margin(10, 10, 10, 10),
        legend.position = "right",
        legend.box = "vertical") +
  ggtitle("STRING PPI colored by module") +
  geom_shadowtext(
    #data = function(d) d[d$name %in% top_hubs, ],
    aes(x = x, y = y, label = label, colour = factor(module)),
    bg.colour = "white",
    bg.r      = 0.1,
    size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) + coord_fixed()

nw_mod

set.seed(27)

# Get layout coordinates
lay <- create_layout(g, layout = "fr")

deg <- igraph::degree(g)

nodes <- lay %>%
  as_tibble() %>%
  mutate(
    label  = as.character(label),
    degree = unname(deg[name])   # name comes from layout, matches vertex name
  )

terms_keep <- c("GO:0006954", "GO:0006898", "GO:0006629") # select the desired terms

go_membership <- enrich_data %>%
  filter(ID %in% terms_keep) %>%
  dplyr::select(ID, geneID) %>%
  mutate(geneID = str_split(geneID, "/")) %>%
  tidyr::unnest(geneID) %>%
  transmute(term = ID, label = geneID) %>%
  distinct()

go_lut <- enrich_data %>%
  dplyr::filter(ID %in% terms_keep) %>%
  dplyr::select(ID, Description) %>%
  dplyr::distinct()

go_labels <- setNames(go_lut$Description, go_lut$ID)
go_labels <- c(go_labels, Other = "Other / not in selected GO terms")

fill_cols <- c("GO:0006954" = "#F15854FF", "GO:0006898" = "#5DA5DAFF", "GO:0006629"= "#FAA43AFF", "Other" = "grey90")

pie_df <- nodes %>%
  dplyr::select(name, x, y, label, degree) %>%
  left_join(go_membership, by = "label") %>%
  mutate(in_term = 1) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = in_term,
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(total = sum(c_across(all_of(terms_keep)))) %>%
  ungroup() %>%
  mutate(across(all_of(terms_keep), ~ ifelse(total > 0, .x / total, NA_real_))) %>%
  mutate(
    Other = ifelse(is.na(total) | total == 0, 1, 0),
    `GO:0006954` = ifelse(is.na(`GO:0006954`), 0, `GO:0006954`),
    `GO:0006898` = ifelse(is.na(`GO:0006898`), 0, `GO:0006898`),
    `GO:0006629` = ifelse(is.na(`GO:0006629`), 0, `GO:0006629`)
  ) %>%
  mutate(r = scales::rescale(degree, to = c(0.2, 0.55)))

set.seed(27)

nw_go <- ggraph(lay) +
  geom_edge_link(aes(alpha = weight), colour = "grey75", show.legend = FALSE) +
  geom_node_point(
    data = pie_df,
    aes(x = x, y = y, size = degree),
    alpha = 0) +
  scale_size_continuous(
    name = "Degree",
    range = c(3, 12),
  ) +
  ggforce::geom_circle(
    data = pie_df,
    aes(x0 = x, y0 = y, r = r),
    inherit.aes = FALSE,
    colour = "black",
    fill = "grey80",
    linewidth = 0.6
  ) +
  scatterpie::geom_scatterpie(
    data = pie_df,
    aes(x = x, y = y, r = r),
    cols = c(terms_keep, "Other"),
    colour = NA, alpha = 0.95
  ) +
  scale_fill_manual(
    values = fill_cols,
    name = "GO term",
    labels = go_labels[c(terms_keep, "Other")]
  ) +
  theme_void() +
  ggtitle("STRING PPI colored by GO terms") +
  geom_shadowtext(
    data = pie_df,
    aes(x = x, y = y, label = label),
    colour = "black", bg.colour = "white", bg.r = 0.12,
    size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  guides(
    fill = guide_legend(order = 1),
    size = guide_legend(
      order = 2,
      override.aes = list(alpha = 1, shape = 21, fill = NA, colour = "black", stroke = 0.4)
    )
  ) +
  theme(legend.position = "right", legend.box = "vertical") + coord_fixed()

nw_go

