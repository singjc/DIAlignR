#' @param mzPntrs (list) a list object containing nested lists of chromatogram data and chromatogram header information
#' @param outfile (char) a character vector of where to write database
#' export
write_mzPntrsdb <- function( mzPntrs, out_file=NULL ){
  if ( is.null(out_file) ){
  out_file <- "./cached_chromatogram_data.mzPntrs"
  }
  
  con <- DBI::dbConnect(RSQLite::SQLite(), out_file)
  
  run_ids <- data.frame(run_id=names(mzPntrs))
  message("Creating run table...")
  head( run_ids )
  DBI::dbCreateTable( conn = con, name = "run", fields = run_ids)
  message("Writing to run table...")
  DBI::dbWriteTable( conn = con, name = "run", value = run_ids, overwrite = T )
  
  chromHead <- data.table::rbindlist(lapply(mzPntrs, `[[`, "chromHead"), idcol = "run_id")
  message("Creating chromHead table...")
  head( chromHead )
  DBI::dbCreateTable( conn = con, name = "chromHead", fields = chromHead )
  DBI::dbWriteTable(conn = con, name = "chromHead", value = chromHead, overwrite = T )
  
  mz <- data.table::rbindlist(lapply(mzPntrs, `[[`, "mz"), idcol = "run_id")
  message("Creating mz table...")
  head( mz )
  DBI::dbCreateTable( conn = con, name = "mz", fields = mz )
  DBI::dbWriteTable( conn = con, name = "mz", value = mz, overwrite = T )
  
  DBI::dbExecute( conn = con, "CREATE INDEX 'chromhead_idx' ON 'chromHead' (
	'run_id',
	'chromatogramId',
	'chromatogramIndex'
);" )
  
  DBI::dbExecute( conn = con, "CREATE INDEX 'mz_idx' ON 'mz' (
	'run_id',
	'CHROMATOGRAM_ID',
	'FRAGMENT_ID',
	'COMPRESSION',
	'DATA_TYPE',
	'DATA'
);"
  )
  
  DBI::dbDisconnect( con )
  
}