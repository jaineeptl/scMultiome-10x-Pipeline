call_macs2_peaks <- function(cfg, obj, group.by = NULL, new_assay = "peaks") {
  if (!isTRUE(cfg$macs2$enabled)) return(obj)

  Seurat::DefaultAssay(obj) <- "ATAC"

  peaks <- Signac::CallPeaks(
    object = obj,
    macs2.path = cfg$macs2$path,
    group.by = group.by
  )

  peaks <- Signac::keepStandardChromosomes(peaks, pruning.mode = "coarse")

  # blacklist differs by genome; Signac includes blacklist_hg38/blacklist_mm10
  blacklist <- if (cfg$species$genome == "hg38") Signac::blacklist_hg38 else Signac::blacklist_mm10
  peaks <- GenomicRanges::subsetByOverlaps(peaks, ranges = blacklist, invert = TRUE)

  mat <- Signac::FeatureMatrix(
    fragments = Signac::Fragments(obj),
    features = peaks,
    cells = colnames(obj)
  )

  # use existing annotation stored in ATAC assay
  ann <- obj[["ATAC"]]@annotation

  obj[[new_assay]] <- Signac::CreateChromatinAssay(
    counts = mat,
    fragments = here::here(cfg$tenx$fragments),
    annotation = ann
  )

  obj
}