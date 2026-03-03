run_multiseq_demux <- function(cfg, obj, meta_col = "MULTIseq_call") {
  if (!isTRUE(cfg$multiseq$enabled)) return(obj)

  stop_if_missing(here::here(cfg$multiseq$lmo_csv), "LMO CSV")
  stop_if_missing(here::here(cfg$multiseq$r1_fastq), "R1 FASTQ")
  stop_if_missing(here::here(cfg$multiseq$r2_fastq), "R2 FASTQ")

  bar.ref <- read.csv(here::here(cfg$multiseq$lmo_csv), header = FALSE)$V1

  # 10x barcodes in Seurat include "-1". MULTIseq often needs first N bases
  cell_ids <- colnames(obj[["RNA"]]@counts)
  cell_ids_trim <- substr(cell_ids, 1, cfg$multiseq$cell_barcode_n)

  rt <- deMULTIplex::MULTIseq.preProcess(
    R1 = here::here(cfg$multiseq$r1_fastq),
    R2 = here::here(cfg$multiseq$r2_fastq),
    cellIDs = cell_ids_trim,
    cell = cfg$multiseq$cell_range,
    umi  = cfg$multiseq$umi_range,
    tag  = cfg$multiseq$tag_range
  )

  bar.table <- deMULTIplex::MULTIseq.align(rt, cell_ids_trim, bar.ref)

  # Save table for inspection
  write_csv(cfg, as.data.frame(bar.table), "multiseq_barcode_table.csv")

  # Simple “top-rank” call (your approach), generalized
  bt <- bar.table
  bt_rank <- apply(bt, 2, rank)
  calls <- apply(bt_rank, 1, function(x) colnames(bt)[which.max(x)])

  # restore "-1" barcode style (Seurat colnames)
  names(calls) <- paste0(names(calls), "-1")

  # attach to Seurat
  obj <- Seurat::AddMetaData(obj, metadata = calls, col.name = meta_col)

  obj
}