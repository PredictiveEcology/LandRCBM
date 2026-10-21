
# Install SpaDES.project
if (tryCatch(packageVersion("SpaDES.project") < "0.1.1", error = function(x) TRUE)){
  install.packages("SpaDES.project", repos = "predictiveecology.r-universe.dev")
}

# Set project path
projectPath <- "~/GitHub/LandRCBM"

# Set up project
out <- SpaDES.project::setupProject(
  Restart = TRUE,
  useGit = "PredictiveEcology",
  overwrite = TRUE,
  paths = list(projectPath = projectPath,
               outputPath  = file.path(projectPath, "outputs", "LandRCBM_SK"),
               modulePath  = file.path(projectPath, "modules"),
               packagePath = file.path(projectPath, "packages"),
               inputPath   = file.path(projectPath, "inputs"),
               cachePath   = file.path(projectPath, "cache")),
  options = options(
    repos = unique(c("predictiveecology.r-universe.dev", getOption("repos"))),
    Require.cloneFrom = Sys.getenv("R_LIBS_USER"),
    spades.moduleCodeChecks = FALSE,
    spades.recoveryMode     = FALSE,
    reproducible.useMemoise = FALSE
  ),
  times = list(start = 2000, end = 2020),
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
      sppEquivCol = "LandR",
      .studyAreaName = "SK",
      dataYear = 2000
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
      standAgesForFitting = c(0, 110),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 200,
      .useCache = "generateData"
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE
    )
  ),
  
  packages = c("reproducible", "googledrive", "sf", "terra", "PredictiveEcology/LandR@development"),
  
  # Align rasterToMatch with SCANFI raster grid
  rasterToMatchTemplate = terra::rast(
    xmin = -2341500, ymin = 5863500, 
    xmax =  3010500, ymax = 9436500,
    res = 30, crs = paste(
      "PROJCRS[\"NAD_1983_Canada_Lambert\",",
      "  BASEGEOGCRS[\"NAD83\",",
      "   DATUM[\"North American Datum 1983\",ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,LENGTHUNIT[\"metre\",1]]],",
      "   PRIMEM[\"Greenwich\",0,ANGLEUNIT[\"degree\",0.0174532925199433]],",
      "   ID[\"EPSG\",4269]],",
      "  CONVERSION[\"Lambert Conic Conformal (2SP)\",",
      "   METHOD[\"Lambert Conic Conformal (2SP)\",ID[\"EPSG\",9802]],",
      "   PARAMETER[\"Latitude of false origin\",0,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8821]],",
      "   PARAMETER[\"Longitude of false origin\",-95,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8822]],",
      "   PARAMETER[\"Latitude of 1st standard parallel\",49,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8823]],",
      "   PARAMETER[\"Latitude of 2nd standard parallel\",77,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8824]],",
      "   PARAMETER[\"Easting at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8826]],",
      "   PARAMETER[\"Northing at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8827]]],",
      "  CS[Cartesian,2],",
      "   AXIS[\"easting\",east,ORDER[1],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]],",
      "   AXIS[\"northing\",north,ORDER[2],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]]]"
    )),
  
  # Study area is SK
  studyArea = {
    reproducible::prepInputs(
      url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
      destinationPath = paths$inputPath
    ) |>
      subset(PRENAME == "Saskatchewan") |>
      sf::st_union() |>
      sf::st_buffer(-120) |>
      sf::st_transform(sf::st_crs(rasterToMatchTemplate)) |>
      sf::st_as_sf()
  },
  studyArea_biomassParam = studyArea,
  rasterToMatch = {
    res <- 240
    terra::rasterize(
      studyArea, 
      terra::crop(rasterToMatchTemplate, studyArea) |>
        terra::aggregate(res / terra::res(rasterToMatchTemplate)[[1]]), 
      touches = FALSE)
  },
  masterRaster = rasterToMatch,
  
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea, dPath = paths$inputPath)
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  }
)
 
out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)

