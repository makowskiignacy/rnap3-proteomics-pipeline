# Refactored protein_analysis_combat.R
# Purpose: Simplify structure, remove repetition, clarify naming, and verify log/linear scale consistency.
# NOTE: Functionality equivalent to original script; outputs saved in same folder structure.

utils::globalVariables(c(
  "PC1","PC2","UMAP1","UMAP2","mean_int","sd_int","CV_percent","p_adj",
  "log2FoldChange","Sig","G1","G2","adj.P.Val","logFC","ProteinId"
))

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(sva)
  library(pheatmap)
  library(RColorBrewer)
  library(ggplot2)
theme_set(theme_bw(base_size = 14))
  library(ggrepel)
  library(umap)
  library(tidyr)
  library(stringr)
  library(readr)
  library(limma)
  library(tibble)
  library(patchwork)
})

# ------------------------------------------------------------
# Utility Functions
# ------------------------------------------------------------
ensure_dir <- function(path) { dir.create(path, recursive = TRUE, showWarnings = FALSE) }

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

safe_source <- function(file) {
  if (file.exists(file)) {
    source(file, local = TRUE)
  } else if (file.exists(file.path("archive", file))) {
    source(file.path("archive", file), local = TRUE)
  } else {
    msg("File", file, "not found (checked ./ and ./archive). Proceeding without it.")
  }
}

# Helper text for transformation descriptions
log_transform_label <- function(original_logged) {
  if (original_logged) "log2 scale (input already log2)" else "log2(x+1) transformed"
}
raw_scale_label <- function(original_logged) {
  if (original_logged) "log2 scale (no original raw provided)" else "raw scale (back-transformed from log2(x+1))"
}

# Heuristic detection if data already appears log-transformed (avoid double log)
# Typical raw proteomics intensities often >> 1e3. If 95th percentile < 100 and max < 1e4 assume it's already on log scale.
is_log_like <- function(x_numeric_vec) {
  x <- as.numeric(x_numeric_vec)
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 10) return(FALSE)
  q95 <- quantile(x, 0.95, na.rm = TRUE)
  mx  <- max(x, na.rm = TRUE)
  # Additional check: many decimals & range not huge
  rng <- diff(range(x, na.rm = TRUE))
  (q95 < 100 && mx < 1e4 && rng < 50)
}

log2_if_needed <- function(mat, label) {
  # Detect on a sample of values
  sample_vals <- as.vector(mat)
  already_log <- is_log_like(sample_vals)
  msg(sprintf("[%s] log-like heuristic: %s", label, ifelse(already_log, "ALREADY LOG", "APPLY log2(x+1)")))
  if (already_log) return(mat) else return(log2(mat + 1))
}

back_transform <- function(log_mat, original_was_logged) {
  if (original_was_logged) {
    # If original data appeared already in log2 scale, do not exponentiate
    return(log_mat)
  } else {
    return((2^log_mat) - 1)
  }
}

quantile_report <- function(mat, label) {
  v <- as.numeric(mat)
  v <- v[is.finite(v)]
  qs <- quantile(v, probs = c(.01,.05,.25,.5,.75,.95,.99), na.rm = TRUE)
  msg(sprintf("%s quantiles: %s", label, paste(names(qs), round(qs,2), collapse = ", ")))
}

# Create z-scores per gene (row-wise), keep only variable genes
zscore_by_row <- function(mat) {
  sds <- apply(mat, 1, sd, na.rm = TRUE)
  keep <- is.finite(sds) & sds > 1e-6
  if (sum(keep) < 2) {
    msg("Z-score heatmap skipped (insufficient variable genes)")
    return(NULL)
  }
  t(scale(t(mat[keep, , drop = FALSE])))
}

save_pheatmap <- function(mat, annotation, annotation_colors, file, title, cluster_cols = TRUE, show_rownames = FALSE) {
  palette <- colorRampPalette(c("#009E73", "white", "#CC79A7"))(100)
  pheatmap(
    mat,
    color = palette,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    main = title,
    show_rownames = show_rownames,
    show_colnames = TRUE,
    cluster_cols = cluster_cols,
    filename = file,
    legend = TRUE,
    fontsize = 12,
    fontsize_row = 12,
    fontsize_col = 12,
    width = 10,
    height = 8
  )
}

run_pca <- function(mat, meta, color_var, shape_var = NULL, title_prefix, out_file, subtitle = NULL, annot_colors = NULL, base_size = 12, width = 4.9, height = 4.0, point_size = 3.5) {
  # Remove rows with any non-finite values
  mat <- mat[apply(mat, 1, function(r) all(is.finite(r))), , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 2) { msg("PCA skipped (insufficient dimensions)"); return(invisible(NULL)) }
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat_f <- mat[vars > 1e-6 & is.finite(vars), , drop = FALSE]
  if (nrow(mat_f) < 2) { msg("PCA skipped (low variance genes)"); return(invisible(NULL)) }
  pcr <- prcomp(t(mat_f), scale. = FALSE)
  pcs <- as.data.frame(pcr$x) %>% rownames_to_column("Sample") %>% left_join(meta, by = c("Sample" = "Sample"))
  pcvar <- summary(pcr)$importance[2, 1:2]
  aes_map <- aes(x = PC1, y = PC2, color = .data[[color_var]])
  if (!is.null(shape_var)) aes_map <- modifyList(aes_map, aes(shape = .data[[shape_var]]))
  g <- ggplot(pcs, aes_map) + geom_point(size = point_size, alpha = 1) +
    labs(title = sprintf("%s (PCA)", title_prefix),
         subtitle = subtitle,
         x = sprintf("PC1: %.1f%%", pcvar[1] * 100),
         y = sprintf("PC2: %.1f%%", pcvar[2] * 100),
         color = color_var, shape = shape_var) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )
  if (!is.null(annot_colors) && color_var %in% names(annot_colors)) {
    g <- g + scale_color_manual(values = annot_colors[[color_var]])
  }
  ggsave(out_file, plot = g, width = width, height = height, dpi = 300)
  g
}

run_umap <- function(mat, meta, color_var, title_prefix, out_file, max_neighbors = 10, subtitle = NULL, annot_colors = NULL, base_size = 12, width = 4.9, height = 4.0, point_size = 3.5) {
  # Remove rows with any non-finite values
  mat <- mat[apply(mat, 1, function(r) all(is.finite(r))), , drop = FALSE]
  if (nrow(mat) < 2 || ncol(mat) < 3) { msg("UMAP skipped (insufficient data)"); return(invisible(NULL)) }
  n_samples <- ncol(mat)
  nn <- min(max_neighbors, n_samples - 1)
  if (nn < 2) { msg("UMAP skipped (neighbors < 2)"); return(invisible(NULL)) }
  # Remove zero / near-zero variance rows
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[vars > 1e-6 & is.finite(vars), , drop = FALSE]
  if (nrow(mat) < 2) { msg("UMAP skipped (low variance genes after filtering)"); return(invisible(NULL)) }
  layout <- umap(t(mat), n_neighbors = nn)$layout
  df <- as.data.frame(layout)
  colnames(df) <- c("UMAP1", "UMAP2")
  df$Sample <- rownames(t(mat))
  df <- left_join(df, meta, by = c("Sample" = "Sample"))
  g <- ggplot(df, aes(UMAP1, UMAP2, color = .data[[color_var]])) +
    geom_point(size = point_size, alpha = 1) + theme_bw(base_size = base_size) +
    labs(title = sprintf("%s (UMAP)", title_prefix), subtitle = subtitle, color = color_var) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    )
  if (!is.null(annot_colors) && color_var %in% names(annot_colors)) {
    g <- g + scale_color_manual(values = annot_colors[[color_var]])
  }
  ggsave(out_file, plot = g, width = width, height = height, dpi = 300)
  g
}

cv_table <- function(long_df, group_vars, intensity_col = "Intensity", min_reps = 2) {
  long_df %>% group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_int = mean(.data[[intensity_col]], na.rm = TRUE),
      sd_int = sd(.data[[intensity_col]], na.rm = TRUE),
      n = sum(!is.na(.data[[intensity_col]])),
      .groups = "drop"
    ) %>%
    filter(n >= min_reps, mean_int > 0, !is.na(sd_int)) %>%
    mutate(CV_percent = (sd_int / mean_int) * 100) %>%
    filter(is.finite(CV_percent))
}

z_sem <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
z_stat <- function(x) mean(x, na.rm = TRUE) / z_sem(x)
z_p <- function(z) (1 - pnorm(abs(z))) * 2

volcano_z <- function(log_mat, comparisons, out_dir, transform_label = NULL) {
  ensure_dir(out_dir)
  get_cols <- function(pattern) grep(paste0("^", pattern, "_\\d+$"), colnames(log_mat), value = TRUE)
  for (nm in names(comparisons)) {
    g1 <- comparisons[[nm]]$group1
    g2 <- comparisons[[nm]]$group2
    c1 <- get_cols(g1)
    c2 <- get_cols(g2)
    if (length(c1) != length(c2) || length(c1) == 0) { msg("Skip", nm, "(unequal replicates)"); next }
    zspace <- log_mat[, c(c1, c2), drop = FALSE] %>% scale() %>% as.data.frame()
    delta <- zspace[, c1] - zspace[, c2]
    zscore <- apply(delta, 1, z_stat)
    pval <- z_p(zscore)
    padj <- p.adjust(pval, method = "BH")
    logfc <- rowMeans(log_mat[, c1], na.rm = TRUE) - rowMeans(log_mat[, c2], na.rm = TRUE)
    res <- data.frame(Gene = rownames(log_mat), log2FoldChange = logfc, z = zscore, p_value = pval, p_adj = padj) %>% arrange(p_adj)
    write.csv(res, file.path(out_dir, paste0("results_", nm, ".csv")), row.names = FALSE)
    res$Sig <- ifelse(res$p_adj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")
    title_txt <- paste("Volcano (z-test)", nm)
    if (!is.null(transform_label)) title_txt <- paste0(title_txt, "\n", transform_label)
    g <- ggplot(res, aes(log2FoldChange, -log10(p_adj), color = Sig)) + geom_point(alpha = 1) +
      scale_color_manual(values = c(Significant = "#CC79A7", `Not Significant` = "#999999")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      labs(title = title_txt, x = "Log2 fold change", y = "-log10(adj p-value)") + theme_bw(base_size = 14)
  ggsave(file.path(out_dir, paste0("volcano_", nm, ".pdf")), plot = g, width = 10, height = 8, dpi = 300)
  }
}

volcano_limma <- function(log_mat, comparisons, out_dir, transform_label = NULL) {
  ensure_dir(out_dir)
  get_cols <- function(pattern) grep(paste0("^", pattern, "_\\d+$"), colnames(log_mat), value = TRUE)
  plot_list <- list()
  for (nm in names(comparisons)) {
    g1 <- comparisons[[nm]]$group1
    g2 <- comparisons[[nm]]$group2
    c1 <- get_cols(g1); c2 <- get_cols(g2)
    comp_cols <- c(c1, c2)
    if (length(c1) == 0 || length(c2) == 0) { msg("Skipping limma (missing cols)", nm); next }
    design <- model.matrix(~ 0 + factor(c(rep("G1", length(c1)), rep("G2", length(c2)))))
    colnames(design) <- c("G1", "G2")
    fit <- lmFit(log_mat[, comp_cols, drop = FALSE], design)
    cont <- makeContrasts(G1 - G2, levels = design)
    fit2 <- eBayes(contrasts.fit(fit, cont))
    res <- topTable(fit2, coef = 1, number = Inf) %>% rownames_to_column("Gene") %>% arrange(adj.P.Val)
    write.csv(res, file.path(out_dir, paste0("results_", nm, ".csv")), row.names = FALSE)
    res$Sig <- ifelse(res$adj.P.Val < 0.05 & abs(res$logFC) > 1, "Significant", "Not Significant")
    title_txt <- paste(nm)
    if (!is.null(transform_label)) title_txt <- paste0(title_txt, "\n", transform_label)
    g <- ggplot(res, aes(logFC, -log10(adj.P.Val), color = Sig)) + geom_point(alpha = 1) +
      scale_color_manual(values = c(Significant = "#CC79A7", `Not Significant` = "#999999")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
      labs(title = title_txt, x = "Log2 fold change", y = "-log10(adj p-value)") +
      theme_bw(base_size = 15.6) +
      theme(
        plot.title = element_text(size = 15.6, face = "bold"),
        axis.title = element_text(size = 15.6, face = "bold"),
        axis.text = element_text(size = 15.6),
        legend.title = element_text(size = 15.6, face = "bold"),
        legend.text = element_text(size = 15.6)
      )
    ggsave(file.path(out_dir, paste0("volcano_", nm, ".pdf")), plot = g, width = 10, height = 8, dpi = 300)
    plot_list[[nm]] <- g
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 2) + plot_layout(guides = "collect")
    ggsave(file.path(out_dir, "combined_volcano.pdf"), plot = combined_plot, width = 13, height = 15, dpi = 300)
  }
}

# ------------------------------------------------------------
# Data Loading & Preparation
# ------------------------------------------------------------
load_id_mapping <- function(mapping_file) {
  id_map <- read_tsv(mapping_file, skip = 1, col_names = c("From", "To"), show_col_types = FALSE)
  msg("ID mappings:", nrow(id_map))
  id_map
}

load_table_s1 <- function(xlsx_file) {
  raw <- read_excel(xlsx_file, skip = 1)
  vals <- raw %>% select(1, 2, 5:16)
  colnames(vals) <- c(
    "ProteinId", "Gene", "W303_s_1", "W303_s_2", "W303_s_3",
    "BY4742_s_1", "BY4742_s_2", "BY4742_s_3",
    "W303_e_1", "W303_e_2", "W303_e_3",
    "BY4742_e_1", "BY4742_e_2", "BY4742_e_3"
  )
  vals %>%
    filter(!is.na(ProteinId), ProteinId != "", !is.na(Gene), Gene != "") %>%
    mutate(across(-c(ProteinId, Gene), ~ readr::parse_number(as.character(.x))))
}

summarise_table_s1 <- function(df) {
  df %>% mutate(
    W303_s = rowMeans(select(., starts_with("W303_s_")), na.rm = TRUE),
    BY4742_s = rowMeans(select(., starts_with("BY4742_s_")), na.rm = TRUE),
    W303_e = rowMeans(select(., starts_with("W303_e_")), na.rm = TRUE),
    BY4742_e = rowMeans(select(., starts_with("BY4742_e_")), na.rm = TRUE)
  ) %>% select(ProteinId, Gene, W303_s, BY4742_s, W303_e, BY4742_e)
}

load_poland_raw <- function(csv_file) {
  read_csv(csv_file, show_col_types = FALSE)
}

summarise_poland_means <- function(df) {
  df %>% select(Accession, `1_1`:`8_3`) %>%
    mutate(
      WT_YPD_Glucose = rowMeans(select(., `1_1`, `1_2`, `1_3`), na.rm = TRUE),
      maf1D_YPD_Glucose = rowMeans(select(., `2_1`, `2_2`, `2_3`), na.rm = TRUE),
      rpc128_YPD_Glucose = rowMeans(select(., `3_1`, `3_2`, `3_3`), na.rm = TRUE),
      WT_YPGly_30C = rowMeans(select(., `4_1`, `4_2`, `4_3`), na.rm = TRUE),
      WT_YPGly_37C = rowMeans(select(., `5_1`, `5_2`, `5_3`), na.rm = TRUE),
      maf1D_YPGly_30C = rowMeans(select(., `6_1`, `6_2`, `6_3`), na.rm = TRUE),
      maf1D_YPGly_37C = rowMeans(select(., `7_1`, `7_2`, `7_3`), na.rm = TRUE),
      rpc128_YPGly_30C = rowMeans(select(., `8_1`, `8_2`, `8_3`), na.rm = TRUE)
    ) %>% select(Accession, WT_YPD_Glucose:rpc128_YPGly_30C)
}

map_accessions <- function(poland_means, id_map) {
  poland_means %>% left_join(id_map, by = c("Accession" = "From")) %>%
    rename(Gene = To) %>% filter(!is.na(Gene), Gene != "") %>% select(-Accession) %>%
    group_by(Gene) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

combine_datasets <- function(tab_s1, poland_unique) {
  comb <- inner_join(summarise_table_s1(tab_s1), poland_unique, by = c("Gene" = "Gene"))
  comb_clean <- comb %>% select(-ProteinId) %>% filter(complete.cases(.))
  list(raw_matrix = {
    m <- as.matrix(comb_clean[,-1]); rownames(m) <- comb_clean$Gene; m
  }, genes = comb_clean$Gene)
}

# Create sample annotation for combined dataset
combined_annotation <- function(poland_unique) {
  pol_cols <- colnames(poland_unique)[colnames(poland_unique) != "Gene"]
  data.frame(
    Sample = c("W303_s", "BY4742_s", "W303_e", "BY4742_e", pol_cols),
    Dataset = c(rep("Table_S1", 4), rep("Poland_Yeast", length(pol_cols))),
    Condition = c("W303_stationary", "BY4742_stationary", "W303_exponential", "BY4742_exponential", pol_cols),
    stringsAsFactors = FALSE
  ) %>% mutate(
    Experimental_Group = case_when(
      str_detect(Condition, "stationary") ~ "Stationary Phase",
      str_detect(Condition, "exponential") ~ "Exponential Phase",
      str_detect(Condition, "YPD_Glucose") ~ "YPD Glucose",
      str_detect(Condition, "YPGly_30C") ~ "YPGlycerol 30°C",
      str_detect(Condition, "YPGly_37C") ~ "YPGlycerol 37°C",
      TRUE ~ "Other"
    )
  )
}

make_annotation_colors <- function(annot) {
  dataset_pal <- c(
    "Table_S1" = "#CC79A7",
    "Poland_Yeast" = "#009E73"
  )
  exp_levels <- unique(annot$Experimental_Group)
  exp_pal_fixed <- c(
    "Stationary Phase" = "#E69F00",
    "Exponential Phase" = "#D55E00",
    "YPD Glucose" = "#009E73",
    "YPGlycerol 30°C" = "#CC79A7",
    "YPGlycerol 30..C" = "#CC79A7",
    "YPGly_30C" = "#CC79A7",
    "YPGlycerol 37°C" = "#0072B2",
    "YPGlycerol 37..C" = "#0072B2",
    "YPGly_37C" = "#0072B2",
    "Other" = "#999999"
  )
  matched_exp <- exp_pal_fixed[names(exp_pal_fixed) %in% exp_levels]
  cond_levels <- sort(unique(annot$Condition))
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7", "#000000")
  cond_pal <- setNames(colorRampPalette(okabe_ito)(length(cond_levels)), cond_levels)
  list(Dataset = dataset_pal, Condition = cond_pal, Experimental_Group = matched_exp)
}

# Batch correction wrapper
run_combat <- function(log_mat, batches) {
  if (length(unique(batches)) < 2) { msg("Single batch – skipping ComBat"); return(log_mat) }
  vars <- apply(log_mat, 1, var, na.rm = TRUE)
  log_mat_f <- log_mat[vars > 1e-6, , drop = FALSE]
  mod <- model.matrix(~ 1, data = data.frame(Batch = batches))
  res <- tryCatch({
    ComBat(dat = log_mat_f, batch = batches, mod = mod, par.prior = TRUE, prior.plots = FALSE)
  }, error = function(e) { msg("ComBat failed:", e$message); log_mat_f })
  res
}

# Polish replicate preparation (returns log matrix with gene rows)
prepare_polish_replicates <- function(poland_raw, id_map) {
  rep <- poland_raw %>% left_join(id_map, by = c("Accession" = "From")) %>% rename(Gene = To) %>%
    filter(!is.na(Gene), Gene != "") %>% select(-Accession) %>%
    group_by(Gene) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  mat <- as.matrix(rep[,-1]); rownames(mat) <- rep$Gene
  # Rename columns to descriptive replicate names (assumes fixed order of 24 columns)
  if (ncol(mat) == 24) {
    colnames(mat) <- c(
      "WT_YPD_Glucose_1", "WT_YPD_Glucose_2", "WT_YPD_Glucose_3",
      "maf1D_YPD_Glucose_1", "maf1D_YPD_Glucose_2", "maf1D_YPD_Glucose_3",
      "rpc128_YPD_Glucose_1", "rpc128_YPD_Glucose_2", "rpc128_YPD_Glucose_3",
      "WT_YPGly_30C_1", "WT_YPGly_30C_2", "WT_YPGly_30C_3",
      "WT_YPGly_37C_1", "WT_YPGly_37C_2", "WT_YPGly_37C_3",
      "maf1D_YPGly_30C_1", "maf1D_YPGly_30C_2", "maf1D_YPGly_30C_3",
      "maf1D_YPGly_37C_1", "maf1D_YPGly_37C_2", "maf1D_YPGly_37C_3",
      "rpc128_YPGly_30C_1", "rpc128_YPGly_30C_2", "rpc128_YPGly_30C_3"
    )
  }
  original_log_like <- is_log_like(mat)
  log_mat <- log2_if_needed(mat, label = "Polish replicate matrix")
  list(log = log_mat, raw_back = back_transform(log_mat, original_log_like), original_logged = original_log_like)
}

polish_annotation <- function(col_names) {
  df <- data.frame(Sample = col_names, stringsAsFactors = FALSE)
  df$Condition <- sub("_(\\d+)$", "", df$Sample)
  df$Replicate <- sub("^.*_(\\d+)$", "\\1", df$Sample)
  df$Strain <- sub("_(YPD|YPGly).*", "", df$Condition)
  df$Experimental_Group <- case_when(
    grepl("YPD_Glucose", df$Condition) ~ "YPD Glucose",
    grepl("YPGly_30C", df$Condition) ~ "YPGlycerol 30°C",
    grepl("YPGly_37C", df$Condition) ~ "YPGlycerol 37°C",
    TRUE ~ "Other"
  )
  df
}

polish_annotation_colors <- function(annot) {
  strains <- sort(unique(annot$Strain))
  okabe_ito <- c("#CC79A7", "#D55E00", "#E69F00", "#56B4E9", "#009E73", "#56B4E9", "#0072B2", "#000000")
  st_cols <- okabe_ito[seq_along(strains)]
  list(
    Strain = setNames(st_cols, strains),
    Experimental_Group = c(
      "YPD Glucose" = "#009E73",
      "YPGlycerol 30°C" = "#CC79A7",
      "YPGlycerol 37°C" = "#0072B2",
      "Other" = "#999999"
    )[unique(annot$Experimental_Group)]
  )
}

# ------------------------------------------------------------
# Main Pipeline
# ------------------------------------------------------------
run_pipeline <- function() {
  msg("=== SESSION INFO ==="); print(sessionInfo())

  # Paths
  protein_file <- "data/Table_S1_wt.xlsx"
  poland_file  <- "data/swathms_mutants.csv"
  mapping_file <- "data/idmapping_2025_03_02.tsv"
  out_root <- "output_plots_combat"
  ensure_dir(out_root)
  dirs <- list(
    combined = file.path(out_root, "combined_analysis"),
    table_s1 = file.path(out_root, "table_s1_analysis"),
    poland = file.path(out_root, "poland_analysis"),
    volcano_z = file.path(out_root, "volcano_plots", "ztest2_pipeline"),
    volcano_limma = file.path(out_root, "volcano_plots", "limma"),
    polish_only = file.path(out_root, "polish_data_analysis")
  ); lapply(dirs, ensure_dir)

  # Load data
  id_map <- load_id_mapping(mapping_file)
  table_s1_vals <- load_table_s1(protein_file)
  poland_raw <- load_poland_raw(poland_file)
  poland_means <- summarise_poland_means(poland_raw)
  poland_unique <- map_accessions(poland_means, id_map)

  # Combine
  combined <- combine_datasets(table_s1_vals, poland_unique)
  combined_raw <- combined$raw_matrix
  quantile_report(combined_raw, "Combined raw means")

  # Log transform with detection
  combined_log <- log2_if_needed(combined_raw, label = "Combined matrix")
  original_combined_was_logged <- is_log_like(combined_raw)
  quantile_report(combined_log, "Combined log scale")
  combined_log_label <- log_transform_label(original_combined_was_logged)
  combined_raw_label <- raw_scale_label(original_combined_was_logged)

  # Annotation
  annot <- combined_annotation(poland_unique)
  annot$Sample <- make.unique(annot$Sample)  # safety
  annot <- annot %>% arrange(Experimental_Group, Condition)
  rownames(annot) <- annot$Sample
  combined_log <- combined_log[, annot$Sample]
  annot_colors <- make_annotation_colors(annot)

  # ComBat
  msg("Running ComBat (combined dataset)...")
  combat_log <- run_combat(combined_log, batches = annot$Dataset)
  combat_log <- combat_log[apply(combat_log, 1, function(r) all(is.finite(r))), , drop = FALSE]
  combat_raw <- back_transform(combat_log, original_combined_was_logged)

  # Diagnostic logs vs raw
  quantile_report(combat_log, "Post-ComBat log")
  quantile_report(combat_raw, "Post-ComBat raw back-transform")
  # Consistency check (only if original not log)
  if (!original_combined_was_logged) {
    if (!isTRUE(all.equal(log2(combat_raw + 1), combat_log, tolerance = 1e-8))) {
      warning("Log / raw mismatch after ComBat (combined)")
    }
  }

  # Heatmaps
  save_pheatmap(combat_raw, annot[, c("Dataset", "Condition", "Experimental_Group")], annot_colors, file = file.path(dirs$combined, "heatmap_combined_data_raw.pdf"), title = sprintf("Combined (%s, ComBat)", combined_raw_label)); save_pheatmap(combat_log, annot[, c("Dataset", "Condition", "Experimental_Group")], annot_colors,
                file = file.path(dirs$combined, "heatmap_combined_data.pdf"),
                title = sprintf("Combined (%s, ComBat)", combined_log_label))
  z_combined <- zscore_by_row(combat_log)
  if (!is.null(z_combined)) {
    save_pheatmap(z_combined, annot[, c("Dataset", "Condition", "Experimental_Group")], annot_colors,
                  file = file.path(dirs$combined, "heatmap_combined_data_zscore.pdf"),
                  title = "Combined (z-score by gene, ComBat)",
                  cluster_cols = TRUE)
  }

  # Boxplots (before & after)
  pre_long <- as.data.frame(combined_log) %>% rownames_to_column("Gene") %>% pivot_longer(-Gene, names_to = "Sample", values_to = "Log2_Int") %>% left_join(annot, by = "Sample")
  post_long <- as.data.frame(combat_log) %>% rownames_to_column("Gene") %>% pivot_longer(-Gene, names_to = "Sample", values_to = "Log2_Int") %>% left_join(annot, by = "Sample")
  bp_fn <- function(df, title, file) {
    g <- ggplot(df, aes(Sample, Log2_Int, fill = Dataset)) + geom_boxplot(outlier.size = 0.5, alpha = 1) +
      scale_fill_manual(values = annot_colors$Dataset) +
      labs(title = title, x = "Sample", y = sprintf("%s intensity", combined_log_label)) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        plot.margin = margin(5, 5, 5, 12)
      )
    ggsave(file, plot = g, width = 4.9, height = 4.0, dpi = 300)
  }
  bp_fn(pre_long, sprintf("%s Before ComBat", combined_log_label), file.path(dirs$combined, "boxplot_intensities_before_combat.pdf"))
  bp_fn(post_long, sprintf("%s After ComBat", combined_log_label), file.path(dirs$combined, "boxplot_intensities_after_combat.pdf"))

  # PCA / UMAP
  run_pca(combat_log, annot, color_var = "Dataset", shape_var = "Experimental_Group",
          title_prefix = "Combined (ComBat)", out_file = file.path(dirs$combined, "PCA_combined_data_by_dataset.pdf"), subtitle = combined_log_label, annot_colors = annot_colors)
  run_pca(combat_log, annot, color_var = "Experimental_Group", shape_var = NULL,
          title_prefix = "Combined by Experimental Group (ComBat)", out_file = file.path(dirs$combined, "PCA_combined_data_by_experimental_group.pdf"), subtitle = combined_log_label, annot_colors = annot_colors)
  run_umap(combat_log, annot, color_var = "Experimental_Group",
           title_prefix = "Combined by Experimental Group (ComBat)", out_file = file.path(dirs$combined, "UMAP_combined_data.pdf"), subtitle = combined_log_label, annot_colors = annot_colors,
           width = 10.0, height = 7.0, base_size = 12)
  run_pca(combat_raw, annot, color_var = "Dataset", shape_var = "Experimental_Group",
          title_prefix = "Combined (ComBat)", out_file = file.path(dirs$combined, "PCA_combined_data_by_dataset_raw.pdf"), subtitle = combined_raw_label, annot_colors = annot_colors)
  run_pca(combat_raw, annot, color_var = "Experimental_Group", shape_var = NULL,
          title_prefix = "Combined by Experimental Group (ComBat)", out_file = file.path(dirs$combined, "PCA_combined_data_by_experimental_group_raw.pdf"), subtitle = combined_raw_label, annot_colors = annot_colors)
  run_umap(combat_raw, annot, color_var = "Experimental_Group",
           title_prefix = "Combined by Experimental Group (ComBat)", out_file = file.path(dirs$combined, "UMAP_combined_data_raw.pdf"), subtitle = combined_raw_label, annot_colors = annot_colors,
           width = 10.0, height = 7.0, base_size = 12)

  # Summary statistics (pre-ComBat combined)
  summary_stats <- as.data.frame(combined_raw) %>% rownames_to_column("Gene") %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Intensity") %>%
    group_by(Sample) %>% summarise(Mean_Intensity = mean(Intensity), Median_Intensity = median(Intensity), SD_Intensity = sd(Intensity), .groups = "drop")
  write.csv(summary_stats, file.path(dirs$combined, "summary_stats_combined_raw.csv"), row.names = FALSE)

  # CV Analysis for Table S1 replicates
  table_long <- table_s1_vals %>% pivot_longer(-c(ProteinId, Gene), names_to = "Sample", values_to = "Intensity") %>%
    separate(Sample, into = c("Strain", "Phase", "Rep"), sep = "_") %>% mutate(Condition = ifelse(Phase == "s", "stationary", ifelse(Phase == "e", "exponential", Phase)))
  table_cv <- cv_table(table_long, group_vars = c("Gene", "Strain", "Condition"))
  # Boxplot with median labels
  med_tbl <- table_cv %>% group_by(Strain, Condition) %>% summarise(Median_CV = median(CV_percent), .groups = "drop") %>% mutate(Group = paste(Strain, Condition, sep = "_"))
  table_cv$Group <- paste(table_cv$Strain, table_cv$Condition, sep = "_")
  table_cv$Group <- factor(table_cv$Group, levels = c("W303_stationary", "BY4742_stationary", "W303_exponential", "BY4742_exponential"))
  group_colors <- c("W303_stationary" = "#CC79A7", "BY4742_stationary" = "#E69F00", "W303_exponential" = "#D55E00", "BY4742_exponential" = "#56B4E9")
  g_cv <- ggplot(table_cv, aes(Group, CV_percent, fill = Group)) + geom_boxplot(alpha = 1) +
    scale_fill_manual(values = group_colors) +
    theme_bw(base_size = 14) + labs(title = sprintf("CV %% by Sample Group"), y = "CV (%)", x = "Condition")
  ggsave(file.path(dirs$table_s1, "cv_analysis_proteins.pdf"), plot = g_cv, width = 10, height = 8, dpi = 300)

  # Poland CV Analysis
  pol_long <- poland_raw %>% left_join(id_map, by = c("Accession" = "From")) %>% rename(Gene = To) %>% filter(!is.na(Gene), Gene != "") %>% select(-Accession) %>%
    group_by(Gene) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-Gene, names_to = "RepCol", values_to = "Intensity") %>%
    mutate(GroupNum = str_split_fixed(RepCol, "_", 2)[,1], Replicate = str_split_fixed(RepCol, "_", 2)[,2]) %>%
    mutate(Condition = recode(GroupNum,
      `1` = "WT_YPD_Glucose", `2` = "maf1D_YPD_Glucose", `3` = "rpc128_YPD_Glucose",
      `4` = "WT_YPGly_30C", `5` = "WT_YPGly_37C", `6` = "maf1D_YPGly_30C",
      `7` = "maf1D_YPGly_37C", `8` = "rpc128_YPGly_30C", .default = "Unknown"))
  pol_cv <- cv_table(pol_long, group_vars = c("Gene", "Condition")) %>% mutate(Experimental_Group = case_when(
    grepl("YPD", Condition) ~ "YPD Glucose", grepl("YPGly_30C", Condition) ~ "YPGlycerol 30°C", grepl("YPGly_37C", Condition) ~ "YPGlycerol 37°C", TRUE ~ "Other"))
  cv_dir <- file.path(dirs$poland, "cv_analysis")
  ensure_dir(cv_dir)
  write.csv(pol_cv, file.path(cv_dir, "poland_cv_analysis_results.csv"), row.names = FALSE)
  pol_cv$Condition <- factor(pol_cv$Condition, levels = c("WT_YPD_Glucose", "maf1D_YPD_Glucose", "rpc128_YPD_Glucose", "WT_YPGly_30C", "maf1D_YPGly_30C", "rpc128_YPGly_30C", "WT_YPGly_37C", "maf1D_YPGly_37C"))
  med_pol <- pol_cv %>% group_by(Condition) %>% summarise(Median_CV = median(CV_percent), .groups = "drop")
  g_pol_cv <- ggplot(pol_cv, aes(Condition, CV_percent, fill = Experimental_Group)) + geom_boxplot(alpha = 1) +
    scale_fill_manual(values = annot_colors$Experimental_Group) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12)
    ) + labs(title = "CV Poland Dataset", y = "CV (%)", x = "Condition", fill = "Growth")
  ggsave(file.path(dirs$poland, "poland_grouped_cv_boxplot.pdf"), plot = g_pol_cv, width = 10.0, height = 6.2, dpi = 300)

  # Polish-only analysis
  polish <- prepare_polish_replicates(poland_raw, id_map)
  pol_log <- polish$log; pol_raw <- polish$raw_back
  polish_original_logged <- polish$original_logged
  polish_log_label <- log_transform_label(polish_original_logged)
  polish_raw_label <- raw_scale_label(polish_original_logged)
  quantile_report(pol_log, "Polish log (heuristic)")
  quantile_report(pol_raw, "Polish raw (heuristic back-transform)")
  pol_annot <- polish_annotation(colnames(pol_log))
  pol_annot_colors <- polish_annotation_colors(pol_annot)
  rownames(pol_annot) <- pol_annot$Sample
  save_pheatmap(pol_raw, pol_annot[, c("Strain", "Experimental_Group")], pol_annot_colors, file.path(dirs$polish_only, "heatmap_polish_data_raw.pdf"), sprintf("SWATH-MS dataset (%s)", polish_raw_label)); save_pheatmap(pol_log, pol_annot[, c("Strain", "Experimental_Group")], pol_annot_colors, file.path(dirs$polish_only, "heatmap_polish_data.pdf"), sprintf("SWATH-MS dataset (%s)", polish_log_label))
  z_polish <- zscore_by_row(pol_log)
  if (!is.null(z_polish)) {
  save_pheatmap(z_polish, pol_annot[, c("Strain", "Experimental_Group")], pol_annot_colors, file.path(dirs$polish_only, "heatmap_polish_data_zscore.pdf"), "SWATH-MS dataset (z-score by gene)", cluster_cols = TRUE)
  }
  run_pca(pol_log, pol_annot, color_var = "Strain", shape_var = "Experimental_Group", title_prefix = "SWATH-MS dataset by Strain", out_file = file.path(dirs$polish_only, "PCA_swathms_data_by_strain_log.pdf"), subtitle = polish_log_label, annot_colors = pol_annot_colors, width = 4.8, height = 4.0, base_size = 12)
  # Added PCA on log-scaled data by Experimental Group
  run_pca(pol_log, pol_annot, color_var = "Experimental_Group", shape_var = NULL, title_prefix = "SWATH-MS dataset by Experimental Group", out_file = file.path(dirs$polish_only, "PCA_swathms_data_by_experimental_group_log.pdf"), subtitle = polish_log_label, annot_colors = pol_annot_colors, width = 4.8, height = 4.0, base_size = 12)
  # Existing PCA on back-transformed (raw) scale
  run_pca(pol_raw, pol_annot, color_var = "Experimental_Group", shape_var = NULL, title_prefix = "SWATH-MS dataset by Group (raw scale)", out_file = file.path(dirs$polish_only, "PCA_polish_data_by_experimental_group_raw.pdf"), subtitle = polish_raw_label, annot_colors = pol_annot_colors)
  run_umap(pol_raw, pol_annot, color_var = "Experimental_Group", title_prefix = "SWATH-MS dataset", out_file = file.path(dirs$polish_only, "UMAP_polish_data_raw.pdf"), subtitle = polish_raw_label, annot_colors = pol_annot_colors)
  # Boxplots
  pol_long_log <- as.data.frame(pol_log) %>% rownames_to_column("Gene") %>% pivot_longer(-Gene, names_to = "Sample", values_to = "Log2_Int") %>% left_join(pol_annot, by = "Sample")
  g_sample <- ggplot(pol_long_log, aes(Sample, Log2_Int, fill = Strain)) + geom_boxplot(alpha = 1) + scale_fill_manual(values = pol_annot_colors$Strain) + theme_bw(base_size = 14) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = sprintf("SWATH-MS dataset %s by Sample", polish_log_label), y = sprintf("%s intensity", polish_log_label))
  ggsave(file.path(dirs$polish_only, "boxplot_polish_data_by_sample.pdf"), plot = g_sample, width = 16, height = 10, dpi = 300)
  g_group <- ggplot(pol_long_log, aes(Experimental_Group, Log2_Int, fill = Experimental_Group)) + geom_boxplot(alpha = 1) + scale_fill_manual(values = pol_annot_colors$Experimental_Group) + theme_bw(base_size = 14) + labs(title = sprintf("SWATH-MS dataset %s by Growth Condition", polish_log_label), y = sprintf("%s intensity", polish_log_label))
  ggsave(file.path(dirs$polish_only, "boxplot_polish_data_by_group.pdf"), plot = g_group, width = 10, height = 8, dpi = 300)
  g_strain <- ggplot(pol_long_log, aes(Strain, Log2_Int, fill = Strain)) + geom_boxplot(alpha = 1) + scale_fill_manual(values = pol_annot_colors$Strain) + theme_bw(base_size = 14) + labs(title = sprintf("SWATH-MS dataset %s by Strain", polish_log_label), y = sprintf("%s intensity", polish_log_label))
  ggsave(file.path(dirs$polish_only, "boxplot_polish_data_by_strain.pdf"), plot = g_strain, width = 8, height = 7, dpi = 300)

  # Volcano plots (z-test & limma) using Polish replicate log matrix
  comparisons <- list(
    "WT_YPD_Glucose_vs_maf1D_YPD_Glucose" = list(group1 = "WT_YPD_Glucose", group2 = "maf1D_YPD_Glucose"),
    "WT_YPD_Glucose_vs_rpc128_YPD_Glucose" = list(group1 = "WT_YPD_Glucose", group2 = "rpc128_YPD_Glucose"),
    "WT_YPGly_30C_vs_maf1D_YPGly_30C" = list(group1 = "WT_YPGly_30C", group2 = "maf1D_YPGly_30C"),
    "WT_YPGly_30C_vs_rpc128_YPGly_30C" = list(group1 = "WT_YPGly_30C", group2 = "rpc128_YPGly_30C"),
    "WT_YPGly_37C_vs_maf1D_YPGly_37C" = list(group1 = "WT_YPGly_37C", group2 = "maf1D_YPGly_37C")
  )
  volcano_z(pol_log, comparisons, dirs$volcano_z, transform_label = polish_log_label)
  volcano_limma(pol_log, comparisons, dirs$volcano_limma, transform_label = polish_log_label)

  msg("=== ANALYSIS COMPLETED (Refactored) ===")
  msg("Diagnostics: If log heatmap values look too low, confirm raw scale quantiles above. If 95th percentile of raw < ~100, data were already log-scale and second log was intentionally skipped.")
}

# Optional: run when script executed directly
if (identical(environment(), globalenv())) {
  try(run_pipeline())
}