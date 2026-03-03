suppressPackageStartupMessages({
  library(here)
  library(yaml)
  library(Seurat)
  library(Signac)
  library(SoupX)
  library(deMULTIplex)
  library(GenomeInfoDb)
  library(ggplot2)
  library(ggpubr)
})

source(here("R/99_utils.R"))
source(here("R/01_config.R"))
source(here("R/02_io_load.R"))
source(here("R/03_multiseq_demux.R"))
source(here("R/04_atac_assay.R"))
source(here("R/05_qc.R"))
source(here("R/06_rna_processing.R"))
source(here("R/07_multimodal_wnn.R"))
source(here("R/08_peak_calling.R"))

cfg <- read_config()
set.seed(cfg$project$seed)

obj <- load_rna_with_optional_soupx(cfg)
obj <- run_multiseq_demux(cfg, obj, meta_col = "multiseq_barcode")

# optional mapping if user provides config/barcode_map.csv
obj <- apply_barcode_map(cfg, obj, meta_col = "multiseq_barcode", out_col = "sample")

obj <- add_atac_assay(cfg, obj)

obj <- rna_qc(cfg, obj)
obj <- atac_qc(cfg, obj)
saveRDS(obj, file.path(cfg$project$out_dir, "checkpoint_post_qc.rds"))

obj <- run_rna_sct(cfg, obj)
obj <- run_rna_umap_cluster(cfg, obj)
saveRDS(obj, file.path(cfg$project$out_dir, "checkpoint_rna.rds"))

# Use ATAC assay directly, or call peaks and use "peaks"
obj <- call_macs2_peaks(cfg, obj, group.by = NULL, new_assay = "peaks")

# If peaks assay exists, run LSI on peaks; else run LSI on ATAC
assay_for_lsi <- if ("peaks" %in% names(obj@assays)) "peaks" else "ATAC"
obj <- run_atac_lsi(cfg, obj, assay = assay_for_lsi)

obj <- run_wnn(cfg, obj)
saveRDS(obj, file.path(cfg$project$out_dir, "multiome_final.rds"))

# generic plots (no hard-coded colors)
p_rna <- DimPlot(obj, reduction = "umap.rna", group.by = "sample", label = FALSE) + ggtitle("RNA UMAP")
save_plot(cfg, "umap_rna.png", p_rna)

p_wnn <- DimPlot(obj, reduction = "wnn.umap", group.by = "sample", label = FALSE) + ggtitle("WNN UMAP")
save_plot(cfg, "umap_wnn.png", p_wnn)