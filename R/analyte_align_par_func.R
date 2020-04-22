analyte_align_par_func <- function( oswdata, mzPntrs, function_param_input ){
  cat("\n---------------------- START --------------------------\n")
  oswdata <- oswdata[[1]]
  mzPntrs <- mzPntrs[[1]]
  function_param_input <- function_param_input[[1]]
  
  analyte <- unique( oswdata$transition_group_id )
  cat( sprintf("analyte: %s\n", analyte) )
  
  # Select reference run based on m-score
  refRunIdx <- getRefRun(oswdata, analyte)
  
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
  
  # Get XIC_group from reference run. if missing, go to next analyte.
  ref <- names(runs)[refRunIdx]
  exps <- setdiff(names(runs), ref)
  chromIndices <- selectChromIndices(oswdata, runname = ref, analyte = analyte, return_index=return_index )
  
  if(is.null(chromIndices)){
    warning("Chromatogram indices for ", analyte, " are missing in ref run ", runs[ref])
    message("Skipping: ", analyte)
    dummy_dt <- data.table::data.table( transition_group_id=NaN, filename=NaN, run_id=NaN, run_type=NaN, alignment_run_id_pair=NaN, alignment_run_file_pair=NaN, RTo=NaN, RT=NaN, leftWidth=NaN, rightWidth=NaN, Intensity=NaN, peak_group_rank=NaN, ms2_m_score=NaN, m_score=NaN )
    dummy_dt$transition_group_id <- analyte
    return( dummy_dt ) #TODO: Return a table maybe?
  } else {
    tictoc::tic()
    XICs.ref <- extractXIC_group(mz = mzPntrs[[ref]]$mz, chromIndices = chromIndices,
                                 XICfilter = function_param_input$XICfilter, SgolayFiltOrd = function_param_input$SgolayFiltOrd,
                                 SgolayFiltLen = function_param_input$SgolayFiltLen)
    ## End timer
    exec_time <- tictoc::toc(quiet = T)
    message(sprintf("Extracting XIC with %s traces for ref run %s: Elapsed Time = %s sec", length(chromIndices), ref, round(exec_time$toc - exec_time$tic, 3) ))
  }
  
  oswdata$run_id_bind <- oswdata$run_id
  
  oswdata %>%
    dplyr::slice( rep(which(run_id==ref & m_score==min(m_score)), each = n() ) ) %>%
    dplyr::mutate( run_id_bind=oswdata$run_id ) %>%
    dplyr::bind_rows( oswdata ) %>%
    dplyr::mutate( run_pair = paste(run_id_bind,ref,sep='_'),
                   run_type = ifelse(run_id==ref, "ref", "eXp")
    ) %>%
    dplyr::group_by( run_pair ) %>%
    tidyr::nest() %>%
    dplyr::filter( run_pair != paste(ref,ref, sep="_") ) %>%
    dplyr::mutate( XICs.ref = list(XICs.ref),
                   mzPntrs = list(mzPntrs),
                   function_param_input = list(function_param_input) ) -> oswdata_runpairs
  
  ## Set-Up for multiple processing
  future::plan( list(future::tweak( future::multiprocess, workers=(future::availableCores()-10) )) )
  
  oswdata_runpairs %>%
    dplyr::mutate( runpair_alignment = furrr::future_pmap( list(data, XICs.ref, mzPntrs, function_param_input), ( ~pairwise_align_par_func(oswdata_runpair_data = data, XICs.ref = XICs.ref, mzPntrs = mzPntrs, function_param_input = function_param_input ) ) ) ) -> tmp
  
  ## Explicitly close multisession workers by switching plan
  future::plan(future::sequential)
  
  analyte_alignment_results <- data.table::rbindlist(tmp$runpair_alignment)
  
  oswdata %>%
    dplyr::filter( run_id==ref & m_score==min(m_score) ) %>%
    dplyr::select( run_id, transition_group_id, filename, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("ms2_m_score"), m_score ) %>%
    dplyr::mutate( run_type = "ref",
                   alignment_run_id_pair = NaN,
                   alignment_run_file_pair = NaN,
                   RTo = RT ) -> ref_data
  
  analyte_alignment_results <- data.table::rbindlist( list(ref_data, analyte_alignment_results), use.names = TRUE, fill = TRUE )
  analyte_alignment_results %>% dplyr::select( transition_group_id, filename, run_id, run_type, alignment_run_id_pair, alignment_run_file_pair, RTo, RT, leftWidth, rightWidth, Intensity, peak_group_rank, contains("ms2_m_score"), m_score ) -> analyte_alignment_results
  
  return( analyte_alignment_results )
}
