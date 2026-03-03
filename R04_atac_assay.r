add_atac_assay <- function(cfg, obj) {
  h5  <- here::here(cfg$tenx$filtered_h5)
  fr  <- here::here(cfg$tenx$fragments)
  pbm <- here::here(cfg$tenx$per_barcode_metrics)

  stop_if_missing(c(h5, fr, pbm), "10x ATAC input")

  counts <- Seurat::Read10X_h5(h5)
  atac_counts <- counts$Peaks %||% counts[["Peaks"]]
  if (is.null(atac_counts)) stop("Could not find ATAC Peaks matrix in filtered_feature_bc_matrix.h5")

  gr <- Signac::StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  keep <- GenomeInfoDb::seqnames(gr) %in% GenomeInfoDb::standardChromosomes(gr)
  atac_counts <- atac_counts[as.vector(keep), ]

  # EnsDb chosen in config
  ensdb_obj <- get(cfg$species$ensdb)
  ann <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_obj)
  GenomeInfoDb::seqlevelsStyle(ann) <- "UCSC"
  GenomeInfoDb::genome(ann) <- cfg$species$genome

  obj[["ATAC"]] <- Signac::CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fr,
    annotation = ann
  )

  # attach per-barcode metrics (peak_region_fragments)
  metrics <- read.csv(pbm)
  rownames(metrics) <- metrics$barcode
  if ("atac_peak_region_fragments" %in% colnames(metrics)) {
    metrics <- metrics[colnames(obj), , drop = FALSE]
    obj$peak_region_fragments <- metrics$atac_peak_region_fragments
  }

  obj
}