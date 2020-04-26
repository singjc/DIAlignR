#' @param db (char) database file with data
#' @param runs (char) vector of run ids to extract data for
#' @param chromIndives (numeric) vector of chromatogram indexes to extract
#' @export
getmzPntrs_on_the_fly <- function( db, runs, chromIndices=NULL ){
  if ( F ){
    db <- "cached_chromatogram_data.mzPntrs"
    runs <- ref
    chromIndices
  }
  
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  message(sprintf( "[DIAlignR::getmzPntrs_on_the_fly] Extracting mzPntrs on the fly for run(s) %s: ", paste(runs, collapse = ", ")), appendLF = FALSE)
  tictoc::tic()
  mzPntrs <- lapply(runs, function(run){
    if( !is.null(chromIndices) ){
      query <- sprintf("select * FROM chromHead where chromHead.run_id = ('%s') and chromHead.chromatogramIndex in (%s)", run, paste(chromIndices, collapse=', '))
    } else {
      query <- sprintf("select * FROM chromHead where chromHead.run_id = ('%s')", run)
    }
    chromHead <- dplyr::collect( dplyr::tbl( con, dbplyr::sql( query )) )
    
    if( !is.null(chromIndices) ){
      query  <- sprintf("select * FROM mz where mz.run_id = ('%s') and mz.CHROMATOGRAM_ID in (%s)", run, paste(chromIndices, collapse=', '))
    } else {
      query <- sprintf("select * FROM mz where mz.run_id = ('%s')", run)
    }
    mz <- dplyr::collect( dplyr::tbl( con, dbplyr::sql( query )) )
    
    list(mz=mz, chromHead=chromHead)
    
  })
  names(mzPntrs) <- runs
  tictoc::toc()
  return( mzPntrs )
}