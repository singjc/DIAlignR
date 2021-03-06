context("Align DIA runs")

test_that("test_alignTargetedRuns",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["XICfilter"]] <- "none"
  params[["globalAlignment"]] <- "loess"
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "sqMass"
  expect_message(
    alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = NULL, applyFun = lapply)
  )
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

  runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
            "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt")
  BiocParallel::register(BiocParallel::MulticoreParam())
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["batchSize"]] <- 10L
  params[["globalAlignment"]] <- "loess"
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "mzML"
  alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = runs, applyFun = BiocParallel::bplapply)
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test2.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  for(i in 4:14){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["XICfilter"]] <- "none"
  params[["context"]] <- "experiment-wide"
  params[["transitionIntensity"]] <- TRUE
  params[["chromFile"]] <- "mzML"
  params[["globalAlignment"]] <- "linear"
  expect_message(
    alignTargetedRuns(dataPath = dataPath,  outFile = "temp", params = params, oswMerged = TRUE,
                      runs = NULL, applyFun = lapply)
  )
  outData <- read.table("temp.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expData <- read.table("test.tsv", stringsAsFactors = FALSE, sep = "\t", header = TRUE)
  expect_identical(dim(outData), dim(expData))
  expect_identical(colnames(outData), colnames(expData))
  expect_identical(outData[["peptide_id"]], expData[["peptide_id"]])
  expect_identical(outData[["precursor"]], expData[["precursor"]])
  expect_identical(outData[["run"]], expData[["run"]])
  x <- sapply(outData[["intensity"]], function(a) sum(as.numeric(strsplit(a, split = ",")[[1]])), USE.NAMES = FALSE)
  expect_equal(x, expData[["intensity"]], tolerance = 1e-04)
  for(i in 6:10){
    expect_equal(outData[[i]], expData[[i]], tolerance = 1e-04)
  }
  file.remove("temp.tsv")

})

test_that("test_getAlignObjs",{
  runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt", "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
  refRun <- "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"
  dataPath <- system.file("extdata", package = "DIAlignR")
  analytes <- c(32L, 898L, 4618L)
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["globalAlignment"]] <- "loess"
  params[["kernelLen"]] <- 13L
  params[["polyOrd"]] <- 4L
  params[["context"]] <- "experiment-wide"
  params[["chromFile"]] <- "mzML"
  expect_warning(
    outData <- getAlignObjs(analytes, runs, dataPath = dataPath, refRun = refRun,
               oswMerged = TRUE, params = params, objType = "light")
    )
  expData <- testAlignObj()
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][[1]], expData, tolerance = 1e-05)
  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["ref"]],
               lapply(XICs[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix), tolerance = 1e-05)
  expect_equal(outData[[2]][["4618"]][["run1_run2"]][["eXp"]],
               lapply(XICs[["hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix), tolerance = 1e-05)
  expData <- data.frame("leftWidth" = 5220.758, "RT" = 5238.35, "rightWidth" = 5261.723)
  expect_equal(as.data.frame(outData[[2]][["4618"]][["run1_run2"]][["peak"]]), expData, tolerance = 1e-05)
  expect_identical(outData[[2]][["32"]], NULL)
  expect_identical(outData[[2]][["898"]], NULL)
})

test_that("test_alignToRef",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["maxPeptideFdr"]] <- 0.05
  params[["context"]] <- "experiment-wide"
  params$kernelLen <- 13L
  params[["globalAlignment"]] <- "linear"
  params[["globalAlignmentFdr"]] <- 0.05
  params[["chromFile"]] <- "mzML"
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  precursors <- getPrecursors(fileInfo, oswMerged= TRUE, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- precursors[precursors$peptide_id %in% c("7040", "9861", "14383"),]

  mzPntrs <- getMZMLpointers(fileInfo)
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  refRuns <- data.frame("peptide_id" = c("7040", "14383", "9861"), "run" = "run1")
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  RSE <- list()
  RSE[["run1_run2"]] <- 38.6594179136227/params$RSEdistFactor
  globalFits <- lapply(globalFits, extractFit, params[["globalAlignment"]])

  data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR, package="DIAlignR")
  data(multipeptide_DIAlignR, package="DIAlignR")

  # Case 1
  df <- multipeptide_DIAlignR[["14383"]]
  df$alignment_rank[2] <- 1L; df$m_score[3] <- 0.06
  XICs.ref <- list()
  XICs.ref[["4618"]] <- lapply(XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt"]][["4618"]], as.matrix)
  outData <- alignToRef(eXp = "run2", ref = "run1", preIdx = 2L, analytes = 4618L, fileInfo, XICs.ref, params, prec2chromIndex,
             mzPntrs, df, globalFits, RSE)
  expect_equal(df[3,], outData[1,])
  expect_equal(outData[2,], data.frame("transition_group_id" = 4618L, feature_id = bit64::NA_integer64_,
                                       RT = 5241.30, intensity = 189.304,  leftWidth = 5224.2, rightWidth = 5265.2, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1, row.names = c(1L)),
               tolerance = 1e-06)

  # Case 2
  df <- multipeptide_DIAlignR[["14383"]]
  outData <- alignToRef(eXp = "run2", ref = "run1", preIdx = 2L, analytes = 4618L, fileInfo, XICs.ref, params, prec2chromIndex,
                        mzPntrs, df, globalFits, RSE)
  df$alignment_rank[3] <- 1L
  expect_equal(outData, df[3,])

  # Case 3
  chromIndices <- prec2chromIndex[["run1"]][c(3,4), "chromatogramIndex"]
  mz <- mzR::openMSfile(file.path(dataPath, "xics","hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML"))
  XICs.ref <- lapply(chromIndices, function(i) extractXIC_group(mz, chromIndices = i))
  names(XICs.ref) <- c("9719", "9720")
  df <- multipeptide_DIAlignR[["9861"]]
  outData <- alignToRef(eXp = "run2", ref = "run1", preIdx = c(3L,4L), analytes = c(9719L, 9720L), fileInfo,
             XICs.ref, params, prec2chromIndex, mzPntrs, df, globalFits, RSE)
  df$alignment_rank[4] <- 1L
  expect_equal(df[4,], outData[1,])
  expect_equal(outData[2,], data.frame("transition_group_id" = 9719L, feature_id = bit64::NA_integer64_,
                                       RT = 2607.05, intensity = 11.80541,  leftWidth = 2591.431, rightWidth = 2625.569, peak_group_rank = NA_integer_,
                                       m_score = NA_real_, run = "run2", alignment_rank = 1, row.names = c(1L)),
               tolerance = 1e-06)

  # Case 4
  expect_message(outData <- alignToRef(eXp = "run2", ref = "run1", preIdx = c(1L), analytes = c(32L), fileInfo,
                        XICs.ref, params, prec2chromIndex, mzPntrs, df, globalFits, RSE))

})

test_that("test_alignIthAnalyte",{
  dataPath <- system.file("extdata", package = "DIAlignR")
  params <- paramsDIAlignR()
  params[["context"]] <- "experiment-wide"
  params$kernelLen <- 13L
  params[["globalAlignmentFdr"]] <- 0.05
  params[["chromFile"]] <- "mzML"
  params[["maxPeptideFdr"]] <- 0.05
  fileInfo <- getRunNames(dataPath, oswMerged = TRUE, params)
  precursors <- getPrecursors(fileInfo, oswMerged= TRUE, params[["runType"]], params[["context"]], params[["maxPeptideFdr"]])
  precursors <- precursors[precursors$peptide_id %in% c("7040", "9861", "14383"),]
  peptideIDs <-  c(7040L, 14383L, 9861L)
  mzPntrs <- getMZMLpointers(fileInfo)
  prec2chromIndex <- getChromatogramIndices(fileInfo, precursors, mzPntrs)
  features <- getFeatures(fileInfo, maxFdrQuery = 0.05)
  multipeptide <- getMultipeptide(precursors, features)
  refRuns <- data.frame("peptide_id" = c("7040", "14383", "9861"), "run" = "run1")
  globalFits <- getGlobalFits(refRuns, features, fileInfo, params[["globalAlignment"]],
                              params[["globalAlignmentFdr"]], params[["globalAlignmentSpan"]])
  RSE <- list()
  RSE[["run1_run2"]] <- RSE[["run1_run0"]] <- 38.6594179136227/params$RSEdistFactor
  globalFits <- lapply(globalFits, extractFit, params[["globalAlignment"]])

  # Case 1
  expect_message(outData <- alignIthAnalyte(rownum = 1, peptideIDs, multipeptide, refRuns, precursors,
                             prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE))
  # Case 2
  outData <- alignIthAnalyte(rownum = 2, peptideIDs, multipeptide, refRuns, precursors,
                             prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE)
  df <- multipeptide[["14383"]]
  df$alignment_rank <- 1L
  df <- df[c(2,1,3),]; row.names(df) <- NULL
  expect_equal(outData, df)
  # Case 3
  outData <- alignIthAnalyte(rownum = 3, peptideIDs, multipeptide, refRuns, precursors,
                             prec2chromIndex, fileInfo, mzPntrs, params, globalFits, RSE)
  expect_equal(outData[6,], data.frame("transition_group_id" = 9719L, feature_id = bit64::NA_integer64_,
                                                     RT = 2607.05, intensity = 11.80541,  leftWidth = 2591.431, rightWidth = 2625.569, peak_group_rank = NA_integer_,
                                                     m_score = NA_real_, run = "run2", alignment_rank = 1, row.names = c(6L)),
                             tolerance = 1e-06)
  expect_identical(outData$alignment_rank , rep(1L, 6))
  expect_equal(outData$intensity[4], 12.92301, tolerance = 1e-06)
})
