% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getRefRun}
\alias{getRefRun}
\title{Fetch the reference run for each peptide}
\usage{
getRefRun(peptideScores, applyFun = lapply)
}
\arguments{
\item{peptideScores}{(list of data-frames) each dataframe has scores of a peptide across all runs.}

\item{applyFun}{(function) value must be either lapply or BiocParallel::bplapply.}
}
\value{
(dataframe) has two columns:
\item{peptide_id}{(integer) a unique id for each peptide.}
\item{run}{(string) run identifier.}
}
\description{
Provides the reference run based on lowest p-value.
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- getRunNames(dataPath = dataPath)
precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_proteomics",
                                context = "experiment-wide", maxPeptideFdr = 0.05)
peptideIDs <- unique(precursorsInfo$peptide_id)
peptidesInfo <- getPeptideScores(fileInfo, peptideIDs)
peptidesInfo <- lapply(peptideIDs, function(pep) dplyr::filter(peptidesInfo, .data$peptide_id == pep))
names(peptidesInfo) <- as.character(peptideIDs)
\dontrun{
getRefRun(peptidesInfo)
}
}
\seealso{
\code{\link{getPeptideScores}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-04-08
}
\keyword{internal}
