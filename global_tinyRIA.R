projectPath <- "~/LandRCBM"
repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
if (!require("SpaDES.project")){
  Require::Install(c("SpaDES.project"), repos = repos, dependencies = TRUE)
}

getRIA <- function(x) {
  x <- sf::st_read(x) 
  ria <- x[x$TSA_NUMBER %in% c('08', '16', '24', '40', '41'),]
  ria <- sf::st_union(ria)|> sf::st_as_sf()
  return(ria)
}

out <- SpaDES.project::setupProject(
  overwrite = FALSE,
  paths = list(projectPath = projectPath,
               outputPath  = file.path(projectPath, "outputs", "tinyRIA"),
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
  times = list(start = 2020, end = 2022),
  modules = c(
    "PredictiveEcology/Biomass_speciesFactorial@development",
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_speciesParameters@development",
    "PredictiveEcology/CBM_defaults@main",
    "PredictiveEcology/Biomass_yieldTables@main",
    "PredictiveEcology/Biomass_core@main",
    "PredictiveEcology/CBM_dataPrep@main",
    "PredictiveEcology/LandRCBM_split3pools@main",
    "PredictiveEcology/CBM_core@main"
  ),
  packages = c("googledrive", 'RCurl', 'XML', "stars", "httr2", "terra"),
  studyArea = {
    reproducible::prepInputs(
      url = "https://drive.google.com/file/d/1LxacDOobTrRUppamkGgVAUFIxNT4iiHU/view?usp=sharing",
      destinationPath = "~/inputs",
      fun = getRIA,
      overwrite = TRUE
    ) |> sf::st_crop(c(xmin = 1000000, xmax = 1200000, ymin = 1100000, ymax = 1300000))
  },
  studyArea_biomassParam = studyArea,
  rasterToMatch = {
    sa <- terra::vect(studyArea)
    targetCRS <- terra::crs(sa)
    rtm <- terra::rast(sa, res = c(250, 250))
    terra::crs(rtm) <- targetCRS
    rtm[] <- 1
    rtm <- terra::mask(rtm, sa)
    rtm
  },
  masterRaster = {
    masterRaster = rasterToMatch
    masterRaster
  },
  rasterToMatch_biomassParam = masterRaster, 
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea,
                                                dPath = "~/inputs")
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  },
  params = list(
    .globals = list(
      dataYear = 2020,
      .plots = c("png"),
      .plotInterval = 10,
      sppEquivCol = 'LandR',
      .studyAreaName = "RIA"
    ),
    Biomass_core = list(.plots = NA),
    Biomass_borealDataPrep = list(
      .studyAreaName = "RIA",
      subsetDataBiomassModel = 50
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 200,
      .plots = "png",
      .useCache = "generateData"
    ),
    Biomass_speciesParameters = list(
      .plots = "png",
      standAgesForFitting = c(0, 110),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE
    )
  )
)

out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)
