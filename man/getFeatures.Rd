% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_osw.R
\name{getFeatures}
\alias{getFeatures}
\title{Get features from all feature files}
\usage{
getFeatures(
  fileInfo,
  maxFdrQuery = 0.05,
  runType = "DIA_proteomics",
  applyFun = lapply
)
}
\arguments{
\item{fileInfo}{(data-frame) output of \code{\link{getRunNames}} function.}

\item{maxFdrQuery}{(numeric) a numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.}

\item{runType}{(char) yhis must be one of the strings "DIA_proteomics", "DIA_Metabolomics".}

\item{applyFun}{(function) value must be either lapply or BiocParallel::bplapply.}
}
\value{
(list of dataframes) each dataframe has following columns:
\item{transition_group_id}{(integer) a unique id for each precursor.}
\item{RT}{(numeric) retention time as in FEATURE.EXP_RT of osw files.}
\item{Intensity}{(numeric) peak intensity as in FEATURE_MS2.AREA_INTENSITY of osw files.}
\item{leftWidth}{(numeric) as in FEATURE.LEFT_WIDTH of osw files.}
\item{rightWidth}{(numeric) as in FEATURE.RIGHT_WIDTH of osw files.}
\item{peak_group_rank}{(integer) rank of each feature associated with transition_group_id.}
\item{m_score}{(numeric) q-value of each feature associated with transition_group_id.}
}
\description{
Get a list of data-frame of OpenSwath features that contains retention time, intensities, boundaries etc.
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- getRunNames(dataPath = dataPath)
features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_proteomics")
dim(features[[2]]) # 938  8
}
\seealso{
\code{\link{getRunNames}, \link{fetchPrecursorsInfo}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2019) + GPL-3
Date: 2019-04-06
}
