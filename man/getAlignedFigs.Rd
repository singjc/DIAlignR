% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualise_chromatograms.R
\name{getAlignedFigs}
\alias{getAlignedFigs}
\title{Plot aligned XICs group for a specific peptide.}
\usage{
getAlignedFigs(
  AlignObj,
  XICs.ref,
  XICs.eXp,
  refPeakLabel,
  annotatePeak = FALSE
)
}
\arguments{
\item{AlignObj}{(S4 object)}

\item{XICs.ref}{(list) List of extracted ion chromatograms (dataframe) from reference run. The dataframe has two columns: first column is for time
and second column indicates intensity.}

\item{XICs.eXp}{(list) List of extracted ion chromatograms (dataframe) from experiment run.The dataframe has two columns: first column is for time
and second column indicates intensity.}

\item{refPeakLabel}{(numeric vector) It contains peak apex, left width and right width.}

\item{annotatePeak}{(logical) TRUE: Peak boundaries and apex will be highlighted.}
}
\value{
A plot to the current device.
}
\description{
AlignObj is the output from getAlignObjs fucntion. This function prepares ggplot objects from AlignObj.
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
runs <- c("hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt",
 "hroest_K120809_Strep10\%PlasmaBiolRepl2_R04_SW_filt")
AlignObjOutput <- getAlignObjs(analytes = 4618L, runs, dataPath = dataPath)
AlignObj <- AlignObjOutput[[2]][["4618"]][[1]][["AlignObj"]]
XICs.ref <- AlignObjOutput[[2]][["4618"]][[1]][["ref"]]
XICs.eXp <- AlignObjOutput[[2]][["4618"]][[1]][["eXp"]]
refPeakLabel <- AlignObjOutput[[2]][["4618"]][[1]][["peak"]]
\dontrun{
getAlignedFigs(AlignObj, XICs.ref, XICs.eXp, refPeakLabel)
}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2019) + GPL-3
Date: 2019-12-13
}
\keyword{internal}
