
# Install SpaDES.project
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)){
  install.packages("SpaDES.project", repos = "predictiveecology.r-universe.dev")
}

# Set project path
projectPath <- "~/LandRCBM"

# Set up project
getRIA <- function(x) {
  x <- sf::st_read(x) 
  ria <- x[x$TSA_NUMBER %in% c('08', '16', '24', '40', '41'),]
  ria <- sf::st_union(ria)|> sf::st_as_sf()
  return(ria)
}

out <- SpaDES.project::setupProject(
  Restart = TRUE,
  useGit = "PredictiveEcology",
  overwrite = TRUE,
  paths = list(projectPath = projectPath,
               outputPath  = file.path(projectPath, "outputs", "tinyRIA"),
               modulePath  = file.path(projectPath, "modules"),
               packagePath = file.path(projectPath, "packages"),
               inputPath   = file.path(projectPath, "inputs"),
               cachePath   = file.path(projectPath, "cache")),
  options = options(
    repos = unique(c("predictiveecology.r-universe.dev", getOption("repos"))),
    Require.cloneFrom = Sys.getenv("R_LIBS_USER"),
    spades.moduleCodeChecks = FALSE,
    spades.recoveryMode = FALSE,
    reproducible.useMemoise = TRUE),
  times = list(start = 2020, end = 2022),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_speciesFactorial@development",
    "PredictiveEcology/Biomass_speciesParameters@development",
    "PredictiveEcology/CBM_defaults@development",
    "PredictiveEcology/Biomass_yieldTables@development",
    "PredictiveEcology/Biomass_core@development",
    "PredictiveEcology/CBM_dataPrep@development",
    "PredictiveEcology/LandRCBM_split3pools@main",
    "PredictiveEcology/CBM_core@development"
  ),
  
  params = list(
    .globals = list(
      .plots = "png",
      .plotInterval = 10,
      .studyAreaName = "RIA",
      dataYear = 2020,
      sppEquivCol = 'LandR'
    ),
    Biomass_borealDataPrep = list(
      subsetDataBiomassModel = 50
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 200,
      .useCache = "generateData"
    ),
    Biomass_speciesParameters = list(
      standAgesForFitting = c(0, 110),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    Biomass_core = list(
      .plots = NA
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE
    )
  ),
  
  packages = c("googledrive", 'RCurl', 'XML', "stars", "httr2", "terra"),
  studyArea = {
    reproducible::prepInputs(
      url = "https://drive.google.com/file/d/1LxacDOobTrRUppamkGgVAUFIxNT4iiHU",
      destinationPath = paths$inputPath,
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
                                                dPath = paths$inputPath)
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  }
)

out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)

