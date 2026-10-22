
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
               outputPath  = file.path(projectPath, "outputs", "LandRCBM_NWT_CanLad"),
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
  times = list(start = 2000, end = 2024),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_speciesFactorial@development",
    "PredictiveEcology/Biomass_speciesParameters@development",
    "PredictiveEcology/CBM_defaults@development",
    "DominiqueCaron/historicalDisturbances@main",
    "PredictiveEcology/Biomass_regeneration@development",
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
      .studyAreaName = "NWT-TaigaPlains",
      dataYear = 2000
    ),
    Biomass_borealDataPrep = list(
      subsetDataBiomassModel = 50,
      adjustAgeAndLongevity = FALSE
    ),
    Biomass_speciesFactorial = list(
      .plots = NULL, #"pdf",
      runExperiment = TRUE,
      factorialSize = "medium"
    ),
    Biomass_speciesParameters = list(
      standAgesForFitting = c(0, 100),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 150,
      .useCache = "generateData"
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE
    ),
    historicalDisturbances = list(
      disturbanceSource = "CanLaD",
      disturbanceTypes = "wildfire"
    )
  ),
  
  packages = c("reproducible", "googledrive", "sf", "terra", "PredictiveEcology/LandR@development"),
  
  # Study area is the Taiga Plains of Northwest Territories
  studyAreaPolygon = {
    nwt <- reproducible::prepInputs(
      url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
      destinationPath = paths$inputPath
    ) |> subset(PRENAME == "Northwest Territories")
    eco <- reproducible::prepInputs(
      url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
      destinationPath = paths$inputPath,
      projectTo = nwt
    ) |> subset(ECOZONE == 4)
    sf::st_intersection(sf::st_union(nwt), sf::st_union(eco)) |> 
      sf::st_as_sf()
  },
  rasterToMatch = {
    # Align rasterToMatch with SCANFI raster grid
    rasterToMatchRes <- 240
    rasterToMatchTemplate <- terra::rast(
      xmin = -2341500, ymin = 5863500,
      xmax =  3010500, ymax = 9436500,
      res = 30, crs = paste(
        "PROJCRS[\"NAD_1983_Canada_Lambert\",",
        "    BASEGEOGCRS[\"NAD83\",",
        "        DATUM[\"North American Datum 1983\",ELLIPSOID[\"GRS 1980\",6378137,298.257222101004,LENGTHUNIT[\"metre\",1]]],",
        "        PRIMEM[\"Greenwich\",0,ANGLEUNIT[\"degree\",0.0174532925199433]],",
        "        ID[\"EPSG\",4269]],",
        "    CONVERSION[\"Lambert Conic Conformal (2SP)\",",
        "        METHOD[\"Lambert Conic Conformal (2SP)\",ID[\"EPSG\",9802]],",
        "        PARAMETER[\"Latitude of false origin\",0,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8821]],",
        "        PARAMETER[\"Longitude of false origin\",-95,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8822]],",
        "        PARAMETER[\"Latitude of 1st standard parallel\",49,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8823]],",
        "        PARAMETER[\"Latitude of 2nd standard parallel\",77,ANGLEUNIT[\"degree\",0.0174532925199433],ID[\"EPSG\",8824]],",
        "        PARAMETER[\"Easting at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8826]],",
        "        PARAMETER[\"Northing at false origin\",0,LENGTHUNIT[\"metre\",1],ID[\"EPSG\",8827]]],",
        "    CS[Cartesian,2],",
        "        AXIS[\"easting\",east,ORDER[1],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]],",
        "        AXIS[\"northing\",north,ORDER[2],LENGTHUNIT[\"metre\",1,ID[\"EPSG\",9001]]]]"
      ))
    studyArea <- sf::st_transform(studyAreaPolygon, sf::st_crs(rasterToMatchTemplate))
    rasterToMatch <- rasterToMatchTemplate |>
      terra::crop(studyArea) |>
      terra::aggregate(rasterToMatchRes / terra::res(rasterToMatchTemplate)[[1]])
    terra::rasterize(studyArea, rasterToMatch)
  },
  studyArea = sf::st_transform(studyAreaPolygon, sf::st_crs(rasterToMatch)),
  studyArea_biomassParam = studyArea,
  masterRaster = rasterToMatch,
  adminLocator = studyArea$PRENAME,
  ecoLocator   = studyArea$ECOZONE,
  
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea, dPath = paths$inputPath)
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  },
  
  outputs = data.frame(
    objectName = c("yieldTablesId", "summaryBySpecies", "rasterToMatch", "disturbanceEvents", "yieldTablesCumulative"),
    file = c("yieldTablesId.rds", "summaryBySpecies.rds", "rasterToMatch.tif", "disturbanceEvents.rds", "yieldTablesCumulative.rds"),
    fun = c("saveRDS", "saveRDS", "writeRaster", "saveRDS", "saveRDS"),
    package = c("base", "base", "terra", "base", "base")
  )
)

out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)

