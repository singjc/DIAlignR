#' @export
pairwise_align_par_func <- function( oswdata_runpair_data, XICs.ref, function_param_input ){
  
  oswdata_runpair_data <- oswdata_runpair_data[[1]]
  XICs.ref <- XICs.ref[[1]]
  # mzPntrs <- mzPntrs[[1]]
  function_param_input <- function_param_input[[1]]
  
  eXp <- subset( oswdata_runpair_data, subset = run_type=="eXp", select = run_id ) %>% unique() %>% as.character()
  ref <- subset( oswdata_runpair_data, subset = run_type=="ref", select = run_id ) %>% unique() %>% as.character()
  analyte <- unique(oswdata_runpair_data$transition_group_id)
  cat(sprintf("Working on analyte: %s | ref: %s | exp: %s", analyte, ref, eXp), file = function_param_input$redirect_output , sep = "\n")
  
  ## Get Overlapping product ms of identifying transitions
  if( function_param_input$identifying ){
    ## Subset list and data into data.table for current analyte
    oswdata_runpair_data %>% dplyr::filter( run_id %in% c(ref, eXp) ) -> oswFile_subset
    
    if ( any(grepl(paste0(".*", function_param_input$filenames$runs[which(rownames(function_param_input$filenames) %in% ref)], ".*"), oswFile_subset$filename)) & any(grepl(paste0(".*", function_param_input$filenames$runs[which(rownames(function_param_input$filenames) %in% eXp)], ".*"), oswFile_subset$filename)) ) {
      
      ## Get overlapping product mz checking for detecting or identifying transition overlap
      procuct_mz_intersect <- Reduce( intersect, lapply(seq(1,dim(oswFile_subset)[1]), function(row_idx) paste( strsplit(oswFile_subset$product_mz, ",")[[row_idx]], strsplit(oswFile_subset$detecting_transitions, ",")[[row_idx]], sep="_" ) ) )
      ## Re-Extract Reference Chrom IDs for only overlapping product mz transitions
      ## TODO need to maybe make this non-repeptitive
      # chromIndices <- selectChromIndices(oswdata, runname = ref, analyte = analyte, product_mz_filter_list=procuct_mz_intersect, return_index=function_param_input$return_index, keep_all_detecting=function_param_input$keep_all_detecting)
      chromIndices <- selectChromIndices(oswdata, runname = ref, analyte = analyte, product_mz_filter_list=procuct_mz_intersect, return_index=function_param_input$return_index, keep_all_detecting=function_param_input$keep_all_detecting)
      tictoc::tic()
      mzPntrs <- getmzPntrs_on_the_fly( db = function_param_input$cached_mzPntrsdb, runs = ref, chromIndices = chromIndices )
      XICs.ref <- extractXIC_group(mz = mzPntrs[[ref]]$mz, chromIndices = chromIndices,
                                   XICfilter = XICfilter, SgolayFiltOrd = SgolayFiltOrd,
                                   SgolayFiltLen = SgolayFiltLen)
      ## End timer
      exec_time <- tictoc::toc(quiet = T)
      cat(sprintf("Re-Extracting XIC with %s traces for ref run %s: Elapsed Time = %s sec", length(chromIndices), ref, round(exec_time$toc - exec_time$tic, 3) ), file = function_param_input$redirect_output , sep = "\n")
    } else {
      ## Set procuct_mz_intersect to NULL
      procuct_mz_intersect <- NULL
      cat( paste0("Chromatogram indices for ", analyte, " are missing in ", function_param_input$runs[eXp]), file = function_param_input$redirect_output , sep = "\n" )
      next
    }
  } else {
    ## Set procuct_mz_intersect to NULL
    procuct_mz_intersect <- NULL
  }
  
  # Get XIC_group from experiment run
  chromIndices <- selectChromIndices(oswdata_runpair_data, runname = eXp, analyte = analyte, product_mz_filter_list=function_param_input$procuct_mz_intersect, return_index=function_param_input$return_index, keep_all_detecting=function_param_input$keep_all_detecting)
  
  if(!is.null(chromIndices)){
    tictoc::tic()
    mzPntrs <- getmzPntrs_on_the_fly( db = function_param_input$cached_mzPntrsdb, runs = eXp, chromIndices = chromIndices )
    cat( sprintf("mzPntrs: %s", mzPntrs), file = function_param_input$redirect_output , sep = "\n" )
    XICs.eXp <- extractXIC_group(mzPntrs[[eXp]]$mz, chromIndices)
    ## End timer
    exec_time <- tictoc::toc(quiet = T)
    cat(sprintf("Extracting XIC with %s traces for eXp run %s: Elapsed Time = %s sec", length(chromIndices), eXp, round(exec_time$toc - exec_time$tic, 3) ), file = function_param_input$redirect_output , sep = "\n")
    # Get the loess fit for hybrid alignment
    loessFits <- list()
    pair <- paste(ref, eXp, sep = "_")
    if(any(pair %in% names(loessFits))){
      Loess.fit <- loessFits[[pair]]
    } else{
      # Loess.fit <- getGlobalAlignment(oswFiles, ref, eXp, maxFdrLoess, spanvalue, fitType = "loess")
      maxFdrLoess_list <- seq(function_param_input$maxFdrLoess, 1, 0.05)
      i <- 1
      Loess.fit <- NULL
      cat( sprintf("Testing maxFdrLoess: %s", maxFdrLoess_i), file = function_param_input$redirect_output , sep = "..." )
      while ( is.null(Loess.fit) & i<length(maxFdrLoess_list) ) {
        maxFdrLoess_i <- maxFdrLoess_list[i]
        Loess.fit <- tryCatch(
          expr = {
            cat( sprintf("%s", maxFdrLoess_i), file = function_param_input$redirect_output , sep = "..." )
            Loess.fit <- getGlobalAlignment(oswdata_runpair_data, ref, eXp, maxFdrLoess_i, function_param_input$spanvalue, fitType = "loess")
          },
          error = function(e){
            # cat(sprintf("ERROR: The following error occured using maxFdrLoess %s: %s", maxFdrLoess_i, e$message), file = function_param_input$redirect_output , sep = "\n")
            Loess.fit <- NULL
          }
        )
        i <- i + 1
        ##TODO Add a stop condition, otherwise loop will for on forever
      }
      tryCatch( expr = {
        cat( sprintf("%s", pair), file = function_param_input$redirect_output , sep = "\n" )
      },
      error = function(e) {
        cat( sprintf("[pairwise_align_par_func] There was an issue with the print statement." ), file = function_param_input$redirect_output , sep = "\n" )
      })
      
      if ( is.null(Loess.fit) ) {
        cat( sprintf("Warn: Was unable to getGlobalAlignment even after permuting different maxFdrLoess thresholds...Skipping...%s", pair), file = function_param_input$redirect_output , sep = "\n" )
        return( NULL ) #TODO change this return to something more representible maybe
      } else {
        cat( sprintf("Successfully got Loess fit using maxFdrLoess: %s", maxFdrLoess_i), file = function_param_input$redirect_output , sep = "\n" ) 
      }
      loessFits[[pair]] <- Loess.fit
      cat( sprintf("Saving Loess Fit Pair for: %s", pair), file = function_param_input$redirect_output , sep = "\n" )
    }
    # Set up constraints for penalizing similarity matrix
    adaptiveRT <- function_param_input$RSEdistFactor*Loess.fit$s
    cat( sprintf("Getting retention time in experiment run mapped to ref run RT"), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("I am HERE..."), file = function_param_input$redirect_output, sep="\n" )
    cat( sprintf("names XICs.ref: %s", names(XICs.ref)), sep="\n", file=function_param_input$redirect_output)
    cat( sprintf("adaptiveRT: %s", adaptiveRT), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("XICs.ref:"), file = function_param_input$redirect_output , sep = "\n" )
    apply(XICs.ref[[1]], 1, function(x){cat(head(x), sep = '\t', file = function_param_input$redirect_output); cat("\n", file=function_param_input$redirect_output)})

    cat( sprintf("dim(XICs.ref): %s\ndim(XICs.eXp): %s\nadaptiveRT: %s\n", paste(unlist(lapply(XICs.ref, function(x) length(x[[1]]))), collapse=", "), paste(unlist(lapply(XICs.eXp, function(x) length(x[[1]]))), collapse=", ")), file = function_param_input$redirect_output , sep = "\n"   )

    # Get retention time in experiment run mapped to reference run retention time.
    eXpRT <- getMappedRT(refRT = function_param_input$refPeak$RT, XICs.ref = XICs.ref, XICs.eXp = XICs.eXp, Loess.fit = Loess.fit, alignType = function_param_input$alignType, adaptiveRT = adaptiveRT, samplingTime = function_param_input$samplingTime, 
                         normalization = function_param_input$normalization, simMeasure = function_param_input$simMeasure, goFactor = function_param_input$goFactor, geFactor = function_param_input$geFactor, cosAngleThresh = function_param_input$cosAngleThresh, 
                         OverlapAlignment = function_param_input$OverlapAlignment, dotProdThresh = function_param_input$dotProdThresh, gapQuantile = function_param_input$gapQuantile, hardConstrain = function_param_input$hardConstrain,
                         samples4gradient = function_param_input$samples4gradient, objType = "light", function_param_input=function_param_input )
    
    tryCatch( expr = {
      cat( sprintf("eXpRT: %s", eXpRT), file = function_param_input$redirect_output , sep = "\n" ) 
    },
    error = function(e) {
      cat( sprintf("[pairwise_align_par_func] There was an issue with the print statement." ), file = function_param_input$redirect_output , sep = "\n" )
    })
    
    
    if ( is.null(eXpRT) ){
      cat( sprintf("[pairwise_align_par_func] eXpRT was null....." ), file = function_param_input$redirect_output , sep = "\n" )
      return( NULL )
    }
    cat( sprintf("Picking Nearest Feature"), file = function_param_input$redirect_output , sep = "\n" )
    eXp_feature <- pickNearestFeature(eXpRT, analyte, oswdata_runpair_data, runname = eXp,
                                      adaptiveRT = adaptiveRT, featureFDR = 0.05)
     
    tryCatch( expr = {
      cat( sprintf("eXp_feature: %s", eXp_feature), file = function_param_input$redirect_output , sep = "\n" )
    },
    error = function(e) {
      cat( sprintf("[pairwise_align_par_func] There was an issue with the print statement." ), file = function_param_input$redirect_output , sep = "\n" )
    })
    cat( sprintf("ref: %s\neXp: %s\nref_run: %s\neXp_run: %s\ndim(XICs.ref): %s\ndim(XICs.eXp): %s\nadaptiveRT: %s\neXpRT: %s\neXp_feature: %s\n", 
                 ref, eXp, 
                 function_param_input$filenames$runs[which(rownames(function_param_input$filenames) %in% ref)], function_param_input$filenames$runs[which(rownames(function_param_input$filenames) %in% eXp)], 
                 paste(unlist(lapply(XICs.ref, function(x) length(x[[1]]))), collapse=", "), paste(unlist(lapply(XICs.eXp, function(x) length(x[[1]]))), collapse=", "), 
                 adaptiveRT, eXpRT, ifelse( is.null(eXp_feature), 'NULL', eXp_feature) ), file = function_param_input$redirect_output , sep = "\n"  )
    
    if(!is.null(eXp_feature)){
      cat("--> Start Writing Results to respective tables", file = function_param_input$redirect_output , sep = "\n")
      # if ( length(eXp_feature) > 1 ) warning( sprintf( "There was more than one feature found!! Taking only first feature.. %s\n", as.character(eXp_feature) ) )
      # lowest m_score index
      if ( any("ms2_m_score" %in% names(eXp_feature)) ){ m_score_filter_name <- "ms2_m_score" } else { m_score_filter_name <- "m_score" }
      use_index <- which( min(eXp_feature[[m_score_filter_name]])==eXp_feature[[m_score_filter_name]] )
      if ( length(use_index)>1 ) { cat("WARN: There were two eXp_features with the same m_score..", file = function_param_input$redirect_output , sep = "\n"); use_index <- 1 }
      # A feature is found. Use this feature for quantification.
      analyte_run_pair_results <- data.table::data.table( transition_group_id = analyte,
                                                          run_id = eXp,
                                                          run_type = "exp",
                                                          alignment_run_id_pair = pair,
                                                          alignment_run_file_pair = paste(unique(oswdata_runpair_data$filename), collapse='_'),
                                                          filename = subset( oswdata_runpair_data, subset = run_type=="eXp", select = filename ) %>% unique %>% as.character(),
                                                          RTo = subset( oswdata_runpair_data, subset = run_type=="eXp" & m_score==min(m_score), select = RT )[[1]][use_index],
                                                          leftWidtho = subset( oswdata_runpair_data, subset = run_type=="eXp" & m_score==min(m_score), select = leftWidth )[[1]][use_index],
                                                          rightWidtho = subset( oswdata_runpair_data, subset = run_type=="eXp" & m_score==min(m_score), select = rightWidth )[[1]][use_index],
                                                          Intensityo = subset( oswdata_runpair_data, subset = run_type=="eXp" & m_score==min(m_score), select = Intensity )[[1]][use_index],
                                                          RT = eXp_feature[["RT"]][use_index],
                                                          leftWidth = eXp_feature[["leftWidth"]][use_index],
                                                          rightWidth = eXp_feature[["rightWidth"]][use_index],
                                                          Intensity = eXp_feature[["Intensity"]][use_index],
                                                          peak_group_rank = eXp_feature[["peak_group_rank"]][use_index],
                                                          d_score = if ( "d_score" %in% names(eXp_feature) ) eXp_feature[["d_score"]][use_index] else "d_score_not_in_table",
                                                          ms2_pep = if ( "ms2_pep" %in% names(eXp_feature) ) eXp_feature[["ms2_pep"]][use_index] else "ms2_pep_not_in_table",
                                                          ms2_m_score = if ( "ms2_m_score" %in% names(eXp_feature) ) eXp_feature[["ms2_m_score"]][use_index] else "see_m_score",
                                                          ipf_pep = if ( "ipf_pep" %in% names(eXp_feature) ) eXp_feature[["ipf_pep"]][use_index] else "no_ipf",
                                                          m_score = eXp_feature[["m_score"]][use_index]
                                                          )
      
      cat("--> Done Writing Results to respective tables\n", file = function_param_input$redirect_output , sep = "\n")
      return( analyte_run_pair_results )
    } else {
      # Feature is not found.
      cat( "WARN: eXp_feature was empty...", file = function_param_input$redirect_output , sep = "\n" )
      return( NULL ) #TODO change this return to something more representible maybe
    }
  } else {
    cat(paste0("Chromatogram indices for ", analyte, " are missing in ", function_param_input$runs[eXp]), file = function_param_input$redirect_output , sep = "\n")
    return( NULL ) #TODO change this return to something more representible maybe
  }
  
}
