## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("transition_group_id", "peak_group_rank", "leftWidth",
                                                        "rightWidth", "RT", "Intensity"))

#' Outputs intensities for each analyte from aligned Targeted-MS runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each analyte.
#' It then align XICs of each analyte to its reference XICs. Best peak, which has lowest m-score, about the aligned retention time is picked for quantification.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @param dataPath (char) Path to mzml and osw directory.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param chrom_ext (char) Extension to search for chromatogram files in data directory. (Default: ".chrom.mzML")
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param maxFdrLoess (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param RSEdistFactor (numeric) This defines how much distance in the unit of rse remains a noBeef zone.
#' @param saveFiles (logical) Must be selected from light, medium and heavy.
#' @param mzPntrs A list of mzRpwiz.
#' @return Two tables of intensity and rention times for every analyte in each run.
#' @seealso \code{\link{getRunNames}, \link{getOswFiles}, \link{getAnalytesName}, \link{getMappedRT}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' intensityTbl <- alignTargetedRuns(dataPath, runs = runs, analytes = c("QFNNTDIVLLEDFQK_3"),
#'  analyteInGroupLabel = FALSE)
#' intensityTbl <- alignTargetedRuns(dataPath, runs = runs, analytes = c("14299_QFNNTDIVLLEDFQK/3"),
#'  analyteInGroupLabel = TRUE)
#' @importFrom dplyr %>%
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
alignTargetedRuns_par <- function(dataPath, alignType = "hybrid", analyteInGroupLabel = FALSE, oswMerged = TRUE,
                                  runs = NULL, analytes = NULL, nameCutPattern = "(.*)(/)(.*)", chrom_ext=".chrom.mzML",
                                  maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 0.01,
                                  spanvalue = 0.1, runType = "DIA_Proteomics",
                                  normalization = "mean", simMeasure = "dotProductMasked",
                                  XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9,
                                  goFactor = 0.125, geFactor = 40,
                                  cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                                  dotProdThresh = 0.96, gapQuantile = 0.5,
                                  hardConstrain = FALSE, samples4gradient = 100,
                                  samplingTime = 3.4,  RSEdistFactor = 3.5, saveFiles = FALSE,
                                  identifying = FALSE, identifying.transitionPEPfilter=0.6, keep_all_detecting=TRUE, mzPntrs = NULL){
  
  if ( F ){
    library(DIAlignR)
    library(mstools)
    library(dplyr)
    library(zoo)
    # dataPath <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/lower_product_mz_threshold/DIAlignR_Analysis/data"
    dataPath <- "/media/roestlab/Data1/User/JustinS/phospho_enriched_u2os/Georges_Results/data"
    dataPath <- "/home/singjust/projects/def-hroest/data/phospho_enriched_U2OS/singjust_results/Georges_Results/data"
    dataPath <- "/media/justincsing/ExtraDrive1/Documents2/Roest_Lab/Github/PTMs_Project/Synth_PhosoPep/Justin_Synth_PhosPep/results/George_lib_repeat2/DIAlignR_Analysis/data/"
    alignType = "hybrid"; analyteInGroupLabel = FALSE; oswMerged = TRUE;
    runs = NULL; analytes = NULL; nameCutPattern = "(.*)(/)(.*)"; chrom_ext=".chrom.sqMass"
    # runs <- c('chludwig_K150309_007b_SW_1_6', 'chludwig_K150309_008_SW_1_4', 'chludwig_K150309_009_SW_1_3', 'chludwig_K150309_010_SW_1_2', 'chludwig_K150309_011_SW_1_1point5', 'chludwig_K150309_012_SW_1_1', 'chludwig_K150309_013_SW_0')
    # runs =  c("chludwig_K150309_013_SW_0", "chludwig_K150309_012_SW_1_1")
    maxFdrQuery = 0.05; maxFdrLoess = 0.01; analyteFDR = 0.01;
    spanvalue = 1; 
    # runType = "DIA_Proteomics";
    runType = "DIA_Proteomics_ipf"
    normalization = "mean"; simMeasure = "dotProductMasked";
    XICfilter = "sgolay"; SgolayFiltOrd = 4; SgolayFiltLen = 9;
    goFactor = 0.125; geFactor = 40;
    cosAngleThresh = 0.3; OverlapAlignment = TRUE;
    dotProdThresh = 0.96; gapQuantile = 0.5;
    hardConstrain = FALSE; samples4gradient = 100;
    samplingTime = 3.4;  RSEdistFactor = 3.5; saveFiles = FALSE;
    mzPntrs = NULL
    identifying=F
    identifying.transitionPEPfilter=0.6
    keep_all_detecting=T
    i=4
    analyteFDR = 1
    analyte<- "ANS(Phospho)SPTTNIDHLK(Label:13C(6)15N(2))_2"
    analyte <- "AGLDNVDAES(Phospho)K(Label:13C(6)15N(2))_2" # Apparently belongs to two peaks
    eXp <- "run11"
    analyte <- "AKNS(Phospho)PEPNEFLR(Label:13C(6)15N(4))_2"
    eXp <- "run5"
    analyte <- "RGS(Phospho)VYHVPLNIVQADAVR(Label:13C(6)15N(4))_3"
    eXp <- "run2"
    analyte <- "DASASS(Phospho)TSTFDAR(Label:13C(6)15N(4))_2"
    ref <- "run12"
    eXp <- "run5"
    analyte <- "RSMS(Phospho)LLGYR(Label:13C(6)15N(4))_2"
    ref <- "run11"
    eXp <- "run0"
    analyte <- "GS(Phospho)VYHVPLNPVQATAVR(Label:13C(6)15N(4))_3"
    ref <- "run11"
    eXp <- "run7"
    maxFdrQuery=1
    maxFdrLoess=0.01
    analyte <- "analyte: .(Acetyl)AAAAAAAGDS(Phospho)DSWDADAFSVEDPVR_3"
    ref <- "run13"
    eXp <- "run5"
  }
  
  
  message(sprintf( "There are %s cores available for multiprocessing, will use %s cores.", future::availableCores(), (future::availableCores()-10) ) )
  
  func_call_start_time <- Sys.time()
  
  ## Save Function input params into a list for parallel input passing
  function_param_input <- list(alignType = alignType, analyteInGroupLabel = analyteInGroupLabel, oswMerged = oswMerged,
                               runs = runs, analytes = analytes, nameCutPattern = nameCutPattern, chrom_ext=chrom_ext,
                               maxFdrQuery = maxFdrQuery, maxFdrLoess = maxFdrLoess, analyteFDR = analyteFDR,
                               spanvalue = spanvalue, runType = runType,
                               normalization = normalization, simMeasure = simMeasure,
                               XICfilter = XICfilter, SgolayFiltOrd = SgolayFiltOrd, SgolayFiltLen = SgolayFiltLen,
                               goFactor = goFactor, geFactor = geFactor,
                               cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                               dotProdThresh = dotProdThresh, gapQuantile = gapQuantile,
                               hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                               samplingTime = samplingTime,  RSEdistFactor = RSEdistFactor, saveFiles = saveFiles,
                               identifying = identifying, identifying.transitionPEPfilter=identifying.transitionPEPfilter, keep_all_detecting=keep_all_detecting) 
  
  # Check if filter length is odd for Savitzky-Golay filter.
  if( (SgolayFiltLen %% 2) != 1){
    return(stop("SgolayFiltLen can only be odd number"))
  }
  
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern, chrom_ext=chrom_ext)
  if(!is.null(runs)){
    filenames <- filenames[filenames$runs %in% runs,]
    missingRun <- setdiff(runs, filenames$runs)
    if(length(missingRun) != 0){
      return(stop(missingRun, " runs are not found."))
    }
  }
  message("Following runs will be aligned:")
  print(filenames[, "runs"], sep = "\n")
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  function_param_input$runs <- runs
  function_param_input$filenames <- filenames
  
  if ( !file.exists(file.path(getwd(), "mzPntrs.rds")) ){

    ## If using mzML files, cache data
    if ( grepl(".*mzML", chrom_ext) ){
      if(is.null(mzPntrs)){
        ######### Collect pointers for each mzML file. #######
        runs <- filenames$runs
        names(runs) <- rownames(filenames)
        # Collect all the pointers for each mzML file.
        message("Collecting metadata from mzML files.")
        # mzPntrs <- getMZMLpointers(dataPath, runs)
        mzPntrs <- getmzPntrs(dataPath, runs)
        message("Metadata is collected from mzML files.")
        return_index <- "chromatogramIndex"
        function_param_input$return_index <- return_index
        ## Save rds object
        saveRDS( mzPntrs, file.path(getwd(), "mzPntrs.rds") )
      }
    } else if ( grepl(".*sqMass", chrom_ext) ){
      if(is.null(mzPntrs)){
        ######### Collect pointers for each mzML file. #######
        runs <- filenames$runs
        names(runs) <- rownames(filenames)
        # Collect all the pointers for each mzML file.
        message("Collecting metadata from sqMass files.")
        # mzPntrs <- getMZMLpointers(dataPath, runs)
        mzPntrs <- getsqMassPntrs(dataPath, runs, nameCutPattern = nameCutPattern, chrom_ext = chrom_ext, .parallel = TRUE)
        message("Metadata is collected from sqMass files.")
        return_index <- "chromatogramIndex"
        function_param_input$return_index <- return_index
        ## Save rds object
        saveRDS( mzPntrs, file.path(getwd(), "mzPntrs.rds") )
      }
    }
  } else {
    message("Reading in a stored mzPntrs.rds object: ", appendLF=FALSE)
    tictoc::tic()
    mzPntrs <- readRDS( file.path(getwd(), "mzPntrs.rds") )
    tictoc::toc()
    return_index <- "chromatogramIndex"
    function_param_input$return_index <- return_index
  }
  
  if ( !file.exists(file.path(getwd(), "oswFiles.rds")) ){
    tictoc::tic()
    oswFiles <- getOswFiles(dataPath, filenames,  maxFdrQuery = maxFdrQuery, analyteFDR = analyteFDR,
                            oswMerged = oswMerged, analytes = NULL, runType = runType, analyteInGroupLabel = analyteInGroupLabel, 
                            identifying = identifying, identifying.transitionPEPfilter=identifying.transitionPEPfilter, mzPntrs = mzPntrs)
    exec_time <- tictoc::toc(quiet = TRUE)
    message( sprintf("[DIAlignR::alignTargetedruns::getOswFiles(#R170)] Extracting OSW results information for %s runs took %s seconds",dim(filenames)[1], round(exec_time$toc - exec_time$tic, 3) ))
    ## Save oswFiles object as rds object
    saveRDS( oswFiles, file.path(getwd(), "oswFiles.rds") )
  } else {
    message("Reading in a stored oswFiles.rds object: ", appendLF=FALSE)
    tictoc::tic()
    oswFiles <- readRDS( file.path(getwd(), "oswFiles.rds") )
    tictoc::toc()
  }
  
  ## Get Reference Analytes
  refAnalytes <- getAnalytesName(oswFiles, analyteFDR, commonAnalytes = FALSE)
  if(!is.null(analytes)){
    analytesFound <- intersect(analytes, refAnalytes)
    analytesNotFound <- setdiff(analytes, analytesFound)
    if(length(analytesNotFound)>0){
      message(paste(analytesNotFound, "not found."))
    }
    refAnalytes <- analytesFound
  }
  
  oswFiles_dt <- data.table::rbindlist(oswFiles, use.names = TRUE, fill = TRUE, idcol = "run_id")
  oswFiles_dt$analytes <- oswFiles_dt$transition_group_id
  
  analyte_row_id_mapping <- data.table::data.table( transition_group_id = unique(oswFiles_dt$transition_group_id) ) %>% dplyr::arrange(transition_group_id) %>% tibble::rowid_to_column(var="analyte_row_id")
  analyte_row_id_mapping$n_analytes <- nrow(analyte_row_id_mapping)
  
  oswFiles_dt %>%
    merge( analyte_row_id_mapping, by = "transition_group_id" ) %>%
    dplyr::group_by( analytes ) %>%
    tidyr::nest() %>%
    dplyr::mutate( mzPntrs = list(mzPntrs),
                   function_param_input = list(function_param_input) ) -> masterTbl
  
  # masterTbl_old <- masterTbl
 masterTbl <- masterTbl_old[c(sample(seq(1, dim(masterTbl_old)[1]), 48)), ] 
  
  
  message("Performing reference-based alignment.")
  start_time <- Sys.time()
  ## Set-Up for multiple processing
  #future::plan( list(future::tweak( future::multicore, workers=10L )) )
  # cl <- future::makeClusterPSOCK(future::availableCores()-30)
  # future::plan(future::cluster, workers = cl)
  
  ## Specify number of workers
  threads <- as.integer(future::availableCores()-4)
  ## Generate worker_id
  worker_id <- rep(1:threads, length.out = nrow(masterTbl))
  ## Add worker ID to data to process
  masterTbl <- bind_cols(tibble(worker_id), masterTbl)
  
  # install.packages("devtools")
  # devtools::install_github("hadley/multidplyr")
  ## Start clusters of n workers
  cluster <- multidplyr::new_cluster( n = threads )
  ## Partition data to send to different workers
  by_worker_id <- masterTbl %>%
    dplyr::group_by( worker_id ) %>%
    multidplyr::partition(., cluster = cluster)


  # Assign libraries
  multidplyr::cluster_library(cluster = cluster, packages = "dplyr")
  multidplyr::cluster_library(cluster = cluster, packages = "DIAlignR")
  # Assign values (use this to load functions or data to each core)
  multidplyr::cluster_copy(cluster = cluster, names = "getMZandChromHead", env = globalenv() )
  
  tictoc::tic()
  by_worker_id %>%
    dplyr::mutate( alignment_results = purrr::pmap( list(data, mzPntrs, function_param_input), ~analyte_align_par_func(oswdata = data, mzPntrs = mzPntrs, function_param_input = function_param_input) ) ) %>%
    collect() %>% # Special collect() function to recombine partitions
    as_tibble() -> tmp
  tictoc::toc()
  rm(by_worker_id)
  gc()
  
  alignment_results <- tryCatch( expr = {
    masterTbl %>%
      dplyr::mutate( alignment_results = purrr::pmap( list(data, mzPntrs, function_param_input), ~analyte_align_par_func(oswdata = data, mzPntrs = mzPntrs, function_param_input = function_param_input) ) ) -> tmp
    
    alignment_results <- data.table::rbindlist(tmp$alignment_results, fill = TRUE)
  }, 
  error = function(e){
    return( data.table::as.data.table(e$message) )
  })
  
  ## Explicitly close multisession workers by switching plan
  ##future::plan(future::sequential)
  # parallel::stopCluster(cl)
  
  # Report the execution time for hybrid alignment step.
  end_time <- Sys.time()
  message("Execution time for alignment = ", end_time - start_time)
  
  func_call_end_time <- Sys.time()
  
  message( sprintf("Function Call started at: %s\nFunction Call ended at: %s\n", func_call_start_time, func_call_end_time) )
  
  ## Cleanup. 
  rm(mzPntrs)
  
  if(saveFiles){
    data.table::fwrite( alignment_results, file = "alignment_results.tsv", sep = "\t" )
    print("Data matrix is available in the current directory")
    return(1)
  } else {
    return(intesityTbl)
  }
}

#' AlignObj for analytes between a pair of runs
#'
#' This function expects osw and mzml directories at dataPath. It first reads osw files and fetches chromatogram indices for each requested analyte.
#' It then align XICs of each analyte to its reference XICs. AlignObj is returned which contains aligned indices and cumulative score along the alignment path.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @param analytes (vector of strings) transition_group_ids for which features are to be extracted. analyteInGroupLabel must be set according the pattern used here.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param refRun (string)
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @param identifying logical value indicating the extraction of identifying transtions. (Default: FALSE)
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param maxFdrLoess (numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in LOESS fit.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param spanvalue (numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param RSEdistFactor (numeric) This defines how much distance in the unit of rse remains a noBeef zone.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @param mzPntrs A list of mzRpwiz.
#' @return A list of AlignObj. Each AlignObj is an S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#' @seealso \code{\link{plotAlignedAnalytes}, \link{getRunNames}, \link{getOswFiles}, \link{getXICs4AlignObj}, \link{getAlignObj}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' AlignObjOutput <- getAlignObjs(analytes = "QFNNTDIVLLEDFQK_3", runs, dataPath = dataPath)
#' plotAlignedAnalytes(AlignObjOutput)
#'
#' @references Gupta S, Ahadi S, Zhou W, Röst H. "DIAlignR Provides Precise Retention Time Alignment Across Distant Runs in DIA and Targeted Proteomics." Mol Cell Proteomics. 2019 Apr;18(4):806-817. doi: https://doi.org/10.1074/mcp.TIR118.001132 Epub 2019 Jan 31.
#'
#' @export
getAlignObjs <- function(analytes, runs, dataPath = ".", alignType = "hybrid",
                         runType = "DIA_Proteomics", refRun = NULL,
                         analyteInGroupLabel = FALSE, identifying = FALSE, oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)",
                         maxFdrQuery = 0.05, maxFdrLoess = 0.01, analyteFDR = 1.00, spanvalue = 0.1,
                         normalization = "mean", simMeasure = "dotProductMasked",
                         XICfilter = "sgolay", SgolayFiltOrd = 4, SgolayFiltLen = 9,
                         goFactor = 0.125, geFactor = 40,
                         cosAngleThresh = 0.3, OverlapAlignment = TRUE,
                         dotProdThresh = 0.96, gapQuantile = 0.5,
                         hardConstrain = FALSE, samples4gradient = 100,
                         samplingTime = 3.4,  RSEdistFactor = 3.5, objType = "light", mzPntrs = NULL){
  
  
  
  if(length(runs) != 2){
    print("For pairwise alignment, two runs are required.")
    return(NULL)
  }
  
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  ##### Get filenames from osw files and check if names are consistent between osw and mzML files. ######
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern, chrom_ext=chrom_ext )
  filenames <- filenames[filenames$runs %in% runs,]
  missingRun <- setdiff(runs, filenames$runs)
  if(length(missingRun) != 0){
    return(stop(missingRun, " runs are not found."))
  }
  
  message("Following runs will be aligned:")
  message(filenames[, "runs"], sep = "\n")
  
  if(is.null(mzPntrs)){
    ######### Collect pointers for each mzML file. #######
    runs <- filenames$runs
    names(runs) <- rownames(filenames)
    # Collect all the pointers for each mzML file.
    message("Collecting metadata from mzML files.")
    mzPntrs <- getMZMLpointers(dataPath, runs)
    message("Metadata is collected from mzML files.")
  }
  
  ######### Get Precursors from the query and respectve chromatogram indices. ######
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = maxFdrQuery, analyteFDR = analyteFDR,
                          oswMerged = oswMerged, analytes = NULL, runType = runType,
                          analyteInGroupLabel = analyteInGroupLabel, identifying = identifying, mzPntrs = mzPntrs)
  
  # Report analytes that are not found
  refAnalytes <- getAnalytesName(oswFiles, analyteFDR, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message(paste(analytesNotFound, "not found."))
  }
  analytes <- analytesFound
  
  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes, XICfilter = XICfilter,
                           SgolayFiltOrd = SgolayFiltOrd, SgolayFiltLen = SgolayFiltLen,
                           mzPntrs = mzPntrs)
  
  ####################### Perfrom alignment ##########################################
  AlignObjs <- vector("list", length(analytes))
  names(AlignObjs) <- analytes
  loessFits <- list()
  message("Perfroming alignment")
  for(analyteIdx in seq_along(analytes)){
    analyte <- analytes[analyteIdx]
    # Select reference run based on m-score
    if(is.null(refRun)){
      refRunIdx <- getRefRun(oswFiles, analyte)
    } else{
      refRunIdx <- which(filenames$runs == refRun)
    }
    
    # Get XIC_group from reference run
    ref <- names(runs)[refRunIdx]
    exps <- setdiff(names(runs), ref)
    XICs.ref <- XICs[[ref]][[analyte]]
    
    # Align experiment run to reference run
    for(eXp in exps){
      # Get XIC_group from experiment run
      XICs.eXp <- XICs[[eXp]][[analyte]]
      if(!is.null(XICs.eXp)){
        # Get the loess fit for hybrid alignment
        pair <- paste(ref, eXp, sep = "_")
        if(any(pair %in% names(loessFits))){
          Loess.fit <- loessFits[[pair]]
        } else{
          Loess.fit <- getGlobalAlignment(oswFiles, ref, eXp, maxFdrLoess, spanvalue, fitType = "loess")
          loessFits[[pair]] <- Loess.fit
        }
        adaptiveRT <-  RSEdistFactor*Loess.fit$s # Residual Standard Error
        # Fetch alignment object between XICs.ref and XICs.eXp
        AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT = adaptiveRT, samplingTime,
                                normalization, simType = simMeasure, goFactor, geFactor,
                                cosAngleThresh, OverlapAlignment,
                                dotProdThresh, gapQuantile, hardConstrain, samples4gradient,
                                objType)
        AlignObjs[[analyte]] <- list()
        # Attach AlignObj for the analyte.
        AlignObjs[[analyte]][[pair]] <- AlignObj
        # Attach intensities of reference XICs.
        AlignObjs[[analyte]][[runs[ref]]] <- XICs.ref
        # Attach intensities of experiment XICs.
        AlignObjs[[analyte]][[runs[eXp]]] <- XICs.eXp
        # Attach peak boundaries to the object.
        ## TODO: Note: Changed peak_group_rank==1 to peak_group_rank==min(peak_group_rank)
        AlignObjs[[analyte]][[paste0(pair, "_pk")]] <- oswFiles[[refRunIdx]] %>%
          dplyr::filter(transition_group_id == analyte & peak_group_rank==min(peak_group_rank)) %>%
          dplyr::select(leftWidth, RT, rightWidth) %>%
          as.vector()
      }
      else {AlignObjs[[analyte]] <- NULL}
    }
  }
  
  ####################### Return AlignedObjs ##########################################
  message("Alignment done. Returning AlignObjs")
  AlignObjs
}
