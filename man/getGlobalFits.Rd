% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_global_fit.R
\name{getGlobalFits}
\alias{getGlobalFits}
\title{Calculates all global alignment needed in refRun}
\usage{
getGlobalFits(
  refRun,
  features,
  fileInfo,
  globalAlignment,
  globalAlignmentFdr,
  globalAlignmentSpan,
  applyFun = lapply
)
}
\arguments{
\item{refRun}{(data-frame) Output of getRefRun function. Must have two columsn : transition_group_id and run.}

\item{features}{(list of data-frames) it is output from getFeatures function.}

\item{fileInfo}{(data-frame) Output of getRunNames function.}

\item{globalAlignment}{(string) Must be from "loess" or "linear".}

\item{globalAlignmentFdr}{(numeric) A numeric value between 0 and 1. Features should have m-score lower than this value for participation in global fit.}

\item{globalAlignmentSpan}{(numeric) Spanvalue for LOESS fit. For targeted proteomics 0.1 could be used.}
}
\value{
(list) Each element is either of class lm or loess.
}
\description{
Calculates all global alignment needed in refRun
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- getRunNames(dataPath, oswMerged = TRUE)
features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
precursors <- getPrecursors(fileInfo, TRUE, "DIA_proteomics", "experiment-wide", 0.01)
precursors <- dplyr::arrange(precursors, .data$peptide_id, .data$transition_group_id)
peptideIDs <- unique(precursors$peptide_id)
peptideScores <- getPeptideScores(fileInfo, peptideIDs, TRUE, "DIA_proteomics", "experiment-wide")
peptideScores <- lapply(peptideIDs, function(pep) dplyr::filter(peptideScores, .data$peptide_id == pep))
names(peptideScores) <- as.character(peptideIDs)
\dontrun{
refRun <- getRefRun(peptideScores)
fits <- getGlobalFits(refRun, features, fileInfo, "linear", 0.05, 0.1)
}
}
\seealso{
\code{\link{getRefRun}, \link{getFeatures}, \link{getGlobalAlignment}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-04-19
}
\keyword{internal}
