#' Outputs AlignObj from an alignment of two XIC-groups
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#'
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param Loess.fit LOESS fit object between reference and experiment run.
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simType (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return A S4 object. Three most-important slots are:
#' \item{indexA_aligned}{(integer) aligned indices of reference run.}
#' \item{indexB_aligned}{(integer) aligned indices of experiment run.}
#' \item{score}{(numeric) cumulative score of alignment.}
#'
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' Loess.fit <- getGlobalAlignment(oswFiles_DIAlignR, ref = "run1", eXp = "run2",
#'  maxFdrGlobal = 0.05, spanvalue = 0.1)
#' AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT = 77.82315,
#'  samplingTime = 3.414, normalization = "mean", simType = "dotProductMasked", goFactor = 0.125,
#'   geFactor = 40, cosAngleThresh = 0.3, OverlapAlignment = TRUE, dotProdThresh = 0.96,
#'   gapQuantile = 0.5, hardConstrain = FALSE, samples4gradient = 100, objType = "light")
#' @export
getAlignObj <- function(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                        normalization, simType, goFactor, geFactor,
                        cosAngleThresh, OverlapAlignment,
                        dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light", function_param_input=NULL){
  if ( !is.null(function_param_input) ){
    cat( sprintf("[getMappedRT -> getAlignObj] Getting noBeef."), file = function_param_input$redirect_output , sep = "\n" )
  }
  # Set up constraints for penalizing similarity matrix
  noBeef <- ceiling(adaptiveRT/samplingTime)
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  B1p <- stats::predict(Loess.fit, tVec.ref[1])
  B2p <- stats::predict(Loess.fit, tVec.ref[length(tVec.ref)])
  if ( !is.null(function_param_input) ){
    cat( sprintf("[getMappedRT -> getAlignObj] Perform dynamic programming for chromatogram alignment"), file = function_param_input$redirect_output , sep = "\n" )
  }
  # Perform dynamic programming for chromatogram alignment
  intensityList.ref <- lapply(XICs.ref, `[[`, 2) # Extracting intensity values
  intensityList.eXp <- lapply(XICs.eXp, `[[`, 2) # Extracting intensity values
  if ( !is.null(function_param_input) ){
    cat( sprintf("[getMappedRT -> getAlignObj] noBeef: %s", noBeef), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] len tVec.ref: %s", length(tVec.ref) ), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] len tVec.eXp: %s", length(tVec.eXp) ), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] B1p: %s", B1p ), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] B2p: %s", B2p ), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] int.ref: %s", intensityList.ref ), file = function_param_input$redirect_output , sep = "\n" )
    cat( sprintf("[getMappedRT -> getAlignObj] int.eXp: %s", intensityList.eXp ), file = function_param_input$redirect_output , sep = "\n" )
  }
  AlignObj <- alignChromatogramsCpp(intensityList.ref, intensityList.eXp,
                                    alignType = "hybrid", tVec.ref, tVec.eXp,
                                    normalization = normalization, simType = simType,
                                    B1p = B1p, B2p = B2p, noBeef = noBeef,
                                    goFactor = goFactor, geFactor = geFactor,
                                    cosAngleThresh = cosAngleThresh, OverlapAlignment = OverlapAlignment,
                                    dotProdThresh = dotProdThresh, gapQuantile = gapQuantile,
                                    hardConstrain = hardConstrain, samples4gradient = samples4gradient,
                                    objType = objType)
  if ( !is.null(function_param_input) ){
    cat( sprintf("[getMappedRT -> getAlignObj] Finished Performing dynamic programming for chromatogram alignment"), file = function_param_input$redirect_output , sep = "\n" )
  }
  AlignObj
}

#' Get mapping of reference RT on experiment run.
#'
#' This function aligns XICs of reference and experiment runs. Using alignment, it maps retention time from refernce run on experiment run.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param refRT Peak's retention-time in reference run.
#' @param XICs.ref List of extracted ion chromatograms from reference run.
#' @param XICs.eXp List of extracted ion chromatograms from experiment run.
#' @param Loess.fit LOESS fit object between reference and experiment run.
#' @param alignType Available alignment methods are "global", "local" and "hybrid".
#' @param adaptiveRT (numeric) Similarity matrix is not penalized within adaptive RT.
#' @param samplingTime (numeric) Time difference between two data-points in each chromatogram. For hybrid and local alignment, samples are assumed to be equally time-spaced.
#' @param normalization (character) Must be selected from "mean", "l2".
#' @param simMeasure (string) Must be selected from dotProduct, cosineAngle,
#' cosine2Angle, dotProductMasked, euclideanDist, covariance and correlation.
#' @param goFactor (numeric) Penalty for introducing first gap in alignment. This value is multiplied by base gap-penalty.
#' @param geFactor (numeric) Penalty for introducing subsequent gaps in alignment. This value is multiplied by base gap-penalty.
#' @param cosAngleThresh (numeric) In simType = dotProductMasked mode, angular similarity should be higher than cosAngleThresh otherwise similarity is forced to zero.
#' @param OverlapAlignment (logical) An input for alignment with free end-gaps. False: Global alignment, True: overlap alignment.
#' @param dotProdThresh (numeric) In simType = dotProductMasked mode, values in similarity matrix higher than dotProdThresh quantile are checked for angular similarity.
#' @param gapQuantile (numeric) Must be between 0 and 1. This is used to calculate base gap-penalty from similarity distribution.
#' @param hardConstrain (logical) If FALSE; indices farther from noBeef distance are filled with distance from linear fit line.
#' @param samples4gradient (numeric) This parameter modulates penalization of masked indices.
#' @param objType (char) Must be selected from light, medium and heavy.
#' @return (numeric)
#' @seealso \code{\link{alignChromatogramsCpp}}
#' @examples
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
#' data(oswFiles_DIAlignR, package="DIAlignR")
#' XICs.ref <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run1"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' XICs.eXp <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run2"]][["14299_QFNNTDIVLLEDFQK/3"]]
#' Loess.fit <- getGlobalAlignment(oswFiles_DIAlignR, ref = "run2", eXp = "run0",
#'  maxFdrGlobal = 0.05, spanvalue = 0.1)
#' adaptiveRT <- 77.82315 #3.5*Loess.fit$s
#' getMappedRT(refRT = 5238.35, XICs.ref, XICs.eXp, Loess.fit, alignType = "hybrid",
#'  adaptiveRT = adaptiveRT, samplingTime = 3.414, normalization = "mean",
#'   simMeasure = "dotProductMasked", goFactor = 0.125, geFactor = 40, cosAngleThresh = 0.3,
#'   OverlapAlignment = TRUE, dotProdThresh = 0.96, gapQuantile = 0.5, hardConstrain = FALSE,
#'   samples4gradient = 100)
#' @export
getMappedRT <- function(refRT, XICs.ref, XICs.eXp, Loess.fit, alignType, adaptiveRT, samplingTime,
                        normalization, simMeasure, goFactor, geFactor, cosAngleThresh,
                        OverlapAlignment, dotProdThresh, gapQuantile, hardConstrain,
                        samples4gradient, objType = "light", function_param_input=NULL ){
  if ( !is.null(function_param_input) ){
    tryCatch( expr = {
      cat( sprintf("[getMappedRT] Getting Align Obj"), file = function_param_input$redirect_output , sep = "\n" )
    },
    error = function(e) {
      cat( sprintf("[getMappedRT] There was an issue with the print statement." ), file = function_param_input$redirect_output , sep = "\n" )
    })
  }
  tryCatch( expr = {
    AlignObj <- getAlignObj(XICs.ref, XICs.eXp, Loess.fit, adaptiveRT, samplingTime,
                            normalization, simType = simMeasure, goFactor, geFactor,
                            cosAngleThresh, OverlapAlignment,
                            dotProdThresh, gapQuantile, hardConstrain, samples4gradient, objType, function_param_input=function_param_input)
  },
  error = function(e){
    if ( !is.null(function_param_input) ){
      
      cat( sprintf("[getMappedRT] Was unabled to get aligned obj:\n%s\nReturning NULL for eXpRT..", e$message), file = function_param_input$redirect_output , sep = "\n" )
      return(NULL)
    }
  })
  
  if ( !is.null(function_param_input) ){
    tryCatch( expr = {
      cat( sprintf("[getMappedRT] Successfully Extracted Align Obj"), file = function_param_input$redirect_output , sep = "\n" )
    },
    error = function(e) {
      cat( sprintf("[getMappedRT] There was an issue with the print statement." ), file = function_param_input$redirect_output , sep = "\n" )
    })
  }
  tVec.ref <- XICs.ref[[1]][["time"]] # Extracting time component
  tVec.eXp <- XICs.eXp[[1]][["time"]] # Extracting time component
  eXpRT <- mappedRTfromAlignObj(refRT, tVec.ref, tVec.eXp, AlignObj)
  
  eXpRT
}

