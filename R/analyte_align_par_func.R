#' @export
analyte_align_par_func <- function( oswdata, function_param_input ){
  
  oswdata <- oswdata[[1]]
  # mzPntrs <- mzPntrs[[1]]
  function_param_input <- function_param_input[[1]]
  
  
  ## Redirect output
  redirect_output <- file( paste("analyte_align_par_func_worker_", unique(oswdata$worker_id_num)[[1]],".txt", sep="" ), open="at" )
  function_param_input$redirect_output <- redirect_output

  
  cat("\n---------------------- START --------------------------\n", file = redirect_output, sep = "\n")
  
  # capture.output(knitr::kable(head(oswdata[,c(1:2,6)])), file=redirect_output, append = T)
  
  analyte <- unique( oswdata$transition_group_id )
  cat( sprintf("Worker %s | Partition Dim %s | analyte (%s of %s): %s\n", 
               unique(oswdata$worker_id_num)[[1]], 
               paste(dim(oswdata), collapse = ", "), 
               unique(oswdata$analyte_row_id), 
               unique(oswdata$n_analytes), 
               analyte), file = redirect_output, sep = "\n" )
  
  ## Outter analyte timer start
  analyte_start_time <- tictoc::tic(quiet = T)
  
  # Select reference run based on m-score
  refRunIdx <- getRefRun(oswdata, analyte)
  ref_run_id <- oswdata$run_id[ refRunIdx ]
  ref_osw_data <- oswdata[refRunIdx, ]
  cat( "refRunIdx\n", refRunIdx, "\n", file=redirect_output, append = T )
  cat( "ref File: ", ref_osw_data$filename[1], "\n", file=redirect_output, append = T )
  oswdata[refRunIdx, ] %>%
    dplyr::group_by( transition_group_id ) %>%
    dplyr::filter(transition_group_id == analyte  & m_score==min(m_score)  ) %>%
    dplyr::ungroup() -> refPeak
  
  ## Do a second pass filter in case there are two features with the same m_score
  ## in that instance, filter on the peakgroup_rank
  if ( dim(refPeak)[1] > 1 ){
    oswdata[refRunIdx, ] %>%
      dplyr::group_by( transition_group_id ) %>%
      dplyr::filter(transition_group_id == analyte  & peak_group_rank==min(peak_group_rank) ) %>%
      dplyr::ungroup() -> refPeak
  }
  
  ## Select useful columns
  refPeak %>%
    dplyr::select(leftWidth, RT, rightWidth, Intensity) -> refPeak
  
  ## Do a third check in case there are still more than one top feature, that happen to have the same RT, peak_group_rank and m_score
  if ( dim(refPeak)[1] > 1 ){
    ## TODO: there may be times when the RT is the same, but the boundaries might be slightly different. Need to account for this somehow
    refPeak <- unique(refPeak)[1,]
  }
  ## Store refPeak for down the line processing comparison
  function_param_input$refPeak <- refPeak
  
  cat(sprintf("Reference Peak at %s with boundaries (%s, %s) and intensity of %s", refPeak$RT, refPeak$leftWidth, refPeak$rightWidth, refPeak$Intensity), file = redirect_output, sep = "\n")
  
  # Get XIC_group from reference run. if missing, go to next analyte.
  ref <- names(function_param_input$runs)[ grepl(ref_run_id, names(function_param_input$runs)) ]
  exps <- setdiff(names(function_param_input$runs), ref)
  chromIndices <- selectChromIndices(oswdata, runname = ref, analyte = analyte, return_index=function_param_input$return_index )
  
  if(is.null(chromIndices)){
    warn_msg <- paste0("Chromatogram indices for ", analyte, " are missing in ref run ", function_param_input$runs[ref])
    warning(warn_msg)
    cat("Skipping: ", analyte, file = redirect_output, sep = "\n")
    dummy_dt <- data.table::data.table( transition_group_id=NaN, filename=NaN, run_id=NaN, run_type=NaN, alignment_run_id_pair=NaN, alignment_run_file_pair=NaN, RTo=NaN, leftWidtho=NaN, rightWidtho=NaN, Intensityo=NaN, RT=NaN, leftWidth=NaN, rightWidth=NaN, Intensity=NaN, peak_group_rank=NaN, ms2_m_score=NaN, m_score=NaN, alignment_log=warn_msg )
    dummy_dt$transition_group_id <- analyte
    return( dummy_dt ) #TODO: Return a table maybe?
  } else {
    tictoc::tic()
    mzPntrs <- getmzPntrs_on_the_fly( db = function_param_input$cached_mzPntrsdb, runs = ref, chromIndices = chromIndices )
    XICs.ref <- extractXIC_group(mz = mzPntrs[[ref]]$mz, chromIndices = chromIndices,
                                 XICfilter = function_param_input$XICfilter, SgolayFiltOrd = function_param_input$SgolayFiltOrd,
                                 SgolayFiltLen = function_param_input$SgolayFiltLen)
    ## End timer
    exec_time <- tictoc::toc(quiet = T)
    cat(sprintf("Extracting XIC with %s traces for ref run %s: Elapsed Time = %s sec", length(chromIndices), ref, round(exec_time$toc - exec_time$tic, 3) ), file = redirect_output, sep = "\n")
  }
  
  oswdata$run_id_bind <- oswdata$run_id
  
  oswdata_runpairs <- tryCatch( expr = {
    
    oswdata %>%
      dplyr::group_by( run_id ) %>%
      dplyr::filter( run_id==ref & m_score==min(m_score) ) %>%
      dplyr::add_count() %>%
      dplyr::filter( ifelse(n>1, ifelse(dplyr::row_number()==1, T, F), T) ) %>%
      dplyr::select( -n ) %>%
      tidyr::uncount( weights = nrow( oswdata ) ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate( run_id_bind = as.matrix(dplyr::select(oswdata, run_id_bind)) )  %>%
      dplyr::bind_rows( oswdata ) %>%
      dplyr::mutate( run_pair = paste(run_id_bind,ref,sep='_'),
                     run_type = ifelse(run_id==ref, "ref", "eXp")
      ) %>%
      dplyr::filter( run_pair != paste(ref,ref, sep="_") ) %>%
      dplyr::arrange( run_pair ) %>%
      dplyr::group_by( run_pair ) %>%
      tidyr::nest()  %>%
      dplyr::mutate( XICs.ref = list(XICs.ref),
                     function_param_input = list(function_param_input) ) -> oswdata_runpairs

  },
  error = function(e){
    cat(sprintf("WARN: There was the following error while setting up oswdata_runpairs for parallel mapping:\n%s", e$message), file = redirect_output, sep = "\n")
    oswdata_runpairs <- NULL
    return(  oswdata_runpairs )
  })
  
  ## Check to make sure oswdata_runpairs is not empty
  if ( nrow(oswdata_runpairs) == 0 ) {
    cat(sprintf("WARN: There was no data for oswdata_runpairs.\n"), file = redirect_output, sep = "\n")
    oswdata_runpairs <- NULL
    return(  oswdata_runpairs )
  }
  
  ## Set-Up for multiple processing
  #future::plan( list(future::tweak( future::multicore, workers=10L )) )
  # cl <- future::makeClusterPSOCK(future::availableCores()-35)
  # future::plan(future::cluster, workers = cl)
  
  tmp_catch <- tryCatch( expr = {
    cat( names(oswdata_runpairs), "\n", file=redirect_output, append = T )
    oswdata_runpairs %>%
      dplyr::mutate( runpair_alignment = purrr::pmap( list(data, XICs.ref, function_param_input), ( ~pairwise_align_par_func(oswdata_runpair_data = data, XICs.ref = XICs.ref, function_param_input = function_param_input ) ) ) ) -> tmp
    cat(sprintf("INFO: Vectorized pairwise alignment was successful!!\n"), file = redirect_output, sep = "\n")
    tmp_catch <- list(data=tmp, alignment_status="Success") 
  },
  error = function(e){
    error_msg <- sprintf("WARN: There was the following error during parallel processing on osw_runpairs using purrr::pmap:\n%s", e$message)
    cat( error_msg, file = redirect_output, sep = "\n" )
    tmp <- NULL
    return( list(data=tmp, alignment_status=error_msg) )
  })
  
  tmp <- tmp_catch$data
  alignment_log <- tmp_catch$alignment_status
  
  
  ## Explicitly close multisession workers by switching plan
  ##future::plan(future::sequential)
  # parallel::stopCluster(cl)
  
  ## Get Reference data
  oswdata %>%
    dplyr::group_by( run_id ) %>%
    dplyr::filter( run_id==ref & m_score==min(m_score) ) %>% 
    dplyr::add_count() %>%
    dplyr::filter( ifelse(n>1, ifelse(dplyr::row_number()==1, T, F), T) ) %>%
    dplyr::select( -n ) %>%
    dplyr::ungroup() %>%
    dplyr::select( run_id, transition_group_id, filename, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("ms2_m_score"), m_score ) %>%
    dplyr::mutate( run_type = "ref",
                   alignment_run_id_pair = NaN,
                   alignment_run_file_pair = NaN,
                   RTo = RT,
                   leftWidtho = leftWidth,
                   rightWidtho = rightWidth,
                   Intensityo = Intensity ) -> ref_data
  
  if ( !is.null(tmp) ){
  analyte_alignment_results <- data.table::rbindlist(tmp$runpair_alignment)
  
  analyte_alignment_results <- data.table::rbindlist( list(ref_data, analyte_alignment_results), use.names = TRUE, fill = TRUE )
  analyte_alignment_results %>% dplyr::select( transition_group_id, filename, run_id, run_type, alignment_run_id_pair, alignment_run_file_pair, RTo, leftWidtho, rightWidtho, Intensityo, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("d_score"), contains("ms2_pep"), contains("ms2_m_score"), contains("ipf_pep"), m_score ) -> analyte_alignment_results
  analyte_alignment_results$alignment_log <- alignment_log
  } else {
    oswdata %>%
      dplyr::group_by( run_id ) %>%
      dplyr::filter( m_score==min(m_score) ) %>%
      dplyr::slice( 1 ) %>%
      dplyr::select( run_id, transition_group_id, filename, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("ms2_m_score"), m_score ) %>%
      dplyr::mutate( run_type = if ( run_id==ref) "ref" else "exp",
                     alignment_run_id_pair = if ( run_id==ref) "NaN" else paste(ref, run_id, sep='_'),
                     alignment_run_file_pair = if ( run_id==ref) "NaN" else paste(ref_data$filename, filename, sep='_'),
                     RTo = RT,
                     leftWidtho = leftWidth,
                     rightWidtho = rightWidth,
                     Intensityo = Intensity ) -> analyte_alignment_results
    analyte_alignment_results$RT <- NaN
    analyte_alignment_results$leftWidth <- NaN
    analyte_alignment_results$rightWidth <- NaN
    analyte_alignment_results$Intensity <- NaN
    analyte_alignment_results %>% dplyr::select( transition_group_id, filename, run_id, run_type, alignment_run_id_pair, alignment_run_file_pair, RTo, leftWidtho, rightWidtho, Intensityo, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("d_score"), contains("ms2_pep"), contains("ms2_m_score"), contains("ipf_pep"), m_score ) -> analyte_alignment_results
    analyte_alignment_results$alignment_log <- alignment_log
    
  }
  # data.table::fwrite( analyte_alignment_results, paste("analyte_alignment_results_align_par_func_worker_", unique(oswdata$worker_id_num)[[1]],".tsv", sep="" ), sep="\t") 
  ## Outer analyte timer end
  analyte_end_time <- tictoc::toc(quiet = T)
  
  cat( sprintf("Done processing analyte (%s of %s): %s -> Elapsed Time = %s sec\n", unique(oswdata$analyte_row_id), unique(oswdata$n_analytes), analyte, round(analyte_end_time$toc - analyte_start_time, 3)), file = redirect_output, sep = "\n" )
  
  return( analyte_alignment_results )
}
