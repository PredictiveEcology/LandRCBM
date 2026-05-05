projectPath <- "~/LandRCBM"
repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
if (!require("SpaDES.project")){
  Require::Install(c("SpaDES.project"), repos = repos, dependencies = TRUE)
}

library("data.table")
library("terra")

out <- SpaDES.project::setupProject(
  overwrite = FALSE,
  paths = list(projectPath = projectPath,
               outputPath  = file.path(projectPath, "outputs", "SK"),
               modulePath  = file.path(projectPath, "modules"),
               packagePath = file.path(projectPath, "packages"),
               inputPath   = file.path(projectPath, "inputs"),
               cachePath   = file.path(projectPath, "cache")),
  options = options(
    repos = c(repos = repos),
    Require.cloneFrom = Sys.getenv("R_LIBS_USER"),
    spades.moduleCodeChecks = FALSE,
    spades.recoveryMode = FALSE,
    reproducible.useMemoise = TRUE),
  times = list(start = 1998, end = 2020),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_speciesFactorial@development",
    "PredictiveEcology/Biomass_speciesParameters@development",
    "PredictiveEcology/CBM_defaults@development",
    "PredictiveEcology/Biomass_yieldTables@main",
    "PredictiveEcology/Biomass_core@development",
    "PredictiveEcology/CBM_dataPrep@development",
    "PredictiveEcology/LandRCBM_split3pools@main",
    "PredictiveEcology/CBM_core@development"
  ),
  packages = c("googledrive", 'RCurl', 'XML', "stars", "httr2", "terra"),
  # Study area is SK
  studyArea = {
    # saskatchewan boundaries
    sk <- reproducible::prepInputs(url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
                                   destinationPath = "inputs")
    sk <- sk[sk$PRENAME == "Saskatchewan",]
    # get ecozones
    eco <- reproducible::prepInputs(url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                                    destinationPath = "inputs", projectTo = sk)
    sa <- reproducible::postProcessTo(eco, cropTo = sk, maskTo = sk) |>
      reproducible::Cache() |>
      sf::st_union() |>
      sf::st_as_sf() |> sf::st_buffer(-125)
    sa
  },
  studyArea_biomassParam = studyArea,
  rasterToMatch = {
    sa <- terra::vect(studyArea)
    rtm <- terra::rast(sa, res = c(250, 250))
    terra::crs(rtm) <- terra::crs(sa)
    sa$id <- 1
    rtm <- terra::rasterize(sa, rtm, field = "id", touches = FALSE)
    rtm
  },
  masterRaster = {
    masterRaster = rasterToMatch
    masterRaster
  },
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea,
                                                dPath = "inputs")
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  },
  params = list(
    .globals = list(
      .plots = c("png"),
      .plotInterval = 10,
      sppEquivCol = 'LandR',
      .studyAreaName = "SK",
      dataYear = 2000
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE
    ),
    Biomass_borealDataPrep = list(
      subsetDataBiomassModel = 50
    ),
    Biomass_speciesFactorial = list(
      .plots = NULL, #"pdf",
      runExperiment = TRUE,
      factorialSize = "medium"
    ),
    Biomass_speciesParameters = list(
      .plots = "png",
      standAgesForFitting = c(0, 110),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 200,
      .plots = "png",
      .useCache = "generateData"
    )
  
  )
  
  out$loadOrder <- unlist(out$modules)
  
  initOut <- SpaDES.core::simInit2(out)
  simOut <- SpaDES.core::spades(initOut)
  
  