#' Get chromatogram data points and chromatogram indexes per chromatogram file
#' @param input A shiny input variable that contains Working Directory Information
#' @param global A list variable containing paths to chromatogram files
#' @return (A list of mzRpwiz)
#' 
#' @importFrom tictoc tic toc
#' 
#' @export
getMZandChromHead <- function( data_table=NULL, chromatogram_file_i=NULL, run=NULL , sql_query=NULL, mzPntrs=NULL ){
  if ( !is.null(data_table) ){
    # print( (data_table) )
    data <- data_table[[1]]
    chromatogram_file_i <- data$chromFiles
    run <- data$run_id
    sql_query <- data$sql_query
  }
  if ( is.null(mzPntrs) ){
    mzPntrs <- list()
  }
  message( sprintf("Processing (%s): %s", run, chromatogram_file_i ) )
  # Establish connection to sqlite database chromatogram file
  conn <- DBI::dbConnect( RSQLite::SQLite(), chromatogram_file_i )
  ## Get table of chromatogram incidces and respective transtion ids
  chromHead <- dplyr::collect( dplyr::tbl(conn, dbplyr::sql(sql_query)) )
  ## Store run id and mz object into master list
  mzPntrs[[run]] <- list()
  ## TODO At somepoint make this store the chrom data for sqmass maybe
  # mzPntrs[[run]]$mz <- mz 
  ## Store file path
  # mzPntrs[[run]]$mz <- chromatogram_file_i
  mzPntrs[[run]]$mz <- dplyr::collect( dplyr::tbl(conn, dbplyr::sql( "SELECT * FROM DATA"  )) )
  mzPntrs[[run]]$chromHead <- chromHead
  ## Append chromHead to sqMass DATA table
  mzPntrs[[run]]$mz <- merge( mzPntrs[[run]]$chromHead , mzPntrs[[run]]$mz, by.x="chromatogramIndex", by.y="CHROMATOGRAM_ID", all=T)
  colnames(mzPntrs[[run]]$mz)[which(grepl("chromatogramIndex", colnames(mzPntrs[[run]]$mz)))] <- "CHROMATOGRAM_ID"
  colnames(mzPntrs[[run]]$mz)[which(grepl("chromatogramId", colnames(mzPntrs[[run]]$mz)))] <- "FRAGMENT_ID"
  
  ## Disconnect form database
  DBI::dbDisconnect(conn) 
  
  if ( !is.null(data) ){
    return( list(mz=mzPntrs[[run]]$mz, chromHead=mzPntrs[[run]]$chromHead) )
  } else {
    return( mzPntrs )
  }
  
}

#' Get a list of data.tables for chromatrogram id mappings
#' @param input A shiny input variable that contains Working Directory Information
#' @param global A list variable containing paths to chromatogram files
#' @return (A list of mzRpwiz)
#' 
#' @importFrom tictoc tic toc
#' 
#' @export
getsqMassPntrs <- function( dataPath, runs, nameCutPattern = "(.*)(/)(.*)", chrom_ext=".chrom.sqMass", .parallel=TRUE  ){
  
  sql_query <- sprintf("SELECT
CHROMATOGRAM.NATIVE_ID AS chromatogramId,
CHROMATOGRAM.ID AS chromatogramIndex
FROM CHROMATOGRAM
")
  
  chromFiles <- list.files(dataPath, pattern = ".sqMass", recursive = T, full.names = T)
  names(chromFiles) <- gsub( chrom_ext, "", basename(chromFiles))
  
  ## Get filenames from osw files and check if names are consistent between osw and mzML files. 
  filenames <- getRunNames( dataPath, oswMerged=TRUE, nameCutPattern = nameCutPattern, chrom_ext = chrom_ext )
  filenames <- filenames[filenames$runs %in% runs,]
  
  tictoc::tic('Pre-Loading mzML Chromatogram Files onto disk')
  if ( .parallel==FALSE ){
    mzPntrs <- list()
    for ( chromatogram_input_index_num in seq(1, length(filenames$runs)) ){
      tryCatch(
        expr = {
          tictoc::tic()
          run <- rownames(filenames)[ chromatogram_input_index_num ]
          current_filename <- filenames$runs[ chromatogram_input_index_num ]
          # message(sprintf("\rCacheing mzML for %s of %s runs", run, length(filenames$runs)))
          ## Get path for current chromatogram file
          chromatogram_file_i <-  chromFiles[ grepl(current_filename, names(chromFiles)) ][[1]]
          mzPntrs <- getMZandChromHead( chromatogram_file_i = chromatogram_file_i, run = run , sql_query = sql_query, mzPntrs = mzPntrs )
          
          ## End timer
          exec_time <- tictoc::toc(quiet = T)
          message(sprintf("\rCacheing sqMass for %s of %s runs: Elapsed Time = %s sec", run, length(filenames$runs), round(exec_time$toc - exec_time$tic, 3) ))
        },
        error = function(e){
          message(sprintf("[getsqMassChromIdMapping] There was an issue cacheing %s, skipping...: %s\n", current_filename, e$message))
        }
      ) # End tryCatch
    }
  } else {
    masterTbl <- merge( data.table::setDT(filenames, keep.rownames = "run_id"), data.table::as.data.table(chromFiles, keep.rownames="runs"), by="runs" )
    
    masterTbl %>%
      dplyr::mutate( sql_query = sql_query ) %>%
      dplyr::select( -filename ) %>%
      dplyr::group_by( runs ) %>%
      tidyr::nest() -> masterTbl
    
    threads <- as.integer(future::availableCores()-2)
    
    # worker_id <- rep(1:threads, length.out = nrow(masterTbl))
    # masterTbl <- bind_cols(tibble(worker_id), masterTbl)
    
    ## Set-Up for multiple processing
    future::plan( list(future::tweak( future::multisession, workers=threads )) )
    
    # # install.packages("devtools")
    # # devtools::install_github("hadley/multidplyr")
    # cluster <- multidplyr::new_cluster( n = threads )
    # 
    # by_worker_id <- masterTbl %>%
    #   dplyr::group_by( worker_id ) %>%
    #   multidplyr::partition(., cluster = cluster)
    # 
    # 
    # # Assign libraries
    # multidplyr::cluster_library(cluster = cluster, packages = "dplyr")
    # # multidplyr::cluster_library("mstools") 
    # # Assign values (use this to load functions or data to each core)
    # multidplyr::cluster_copy(cluster = cluster, names = "getMZandChromHead", env = globalenv() )
    # 
    
    
    # tictoc::tic()
    # by_worker_id %>%
    # dplyr::mutate( mzPntrs = purrr::pmap( list(data), ~getMZandChromHead( data_table =data ), .options = furrr::future_options(globals = "getMZandChromHead", scheduling = 2 )) ) %>%
    #   collect() %>% # Special collect() function to recombine partitions
    #   as_tibble() -> tmp
    # tictoc::toc()
    # rm(by_worker_id)
    # gc()
    tictoc::tic()
    masterTbl %>%
      dplyr::mutate( mzPntrs = furrr::future_pmap( list(data), ~getMZandChromHead( data_table = data ), .options = furrr::future_options(globals = "getMZandChromHead", scheduling = 2 )) ) -> tmp
    tictoc::toc()
    
    ## Explicitly close multisession workers by switching plan
    future::plan(future::sequential)
    
    
    mzPntrs <- tmp$mzPntrs
    names(mzPntrs) <- tidyr::unnest(masterTbl, cols="data") %>% dplyr::ungroup() %>% dplyr::select( run_id ) %>% as.matrix() %>% as.character()
    
    
  }
  
  tictoc::toc()
  
  return( mzPntrs )
}