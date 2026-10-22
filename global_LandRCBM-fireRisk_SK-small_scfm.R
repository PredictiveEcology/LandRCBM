
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
               outputPath  = file.path(projectPath, "outputs", "LandRCBM-fireRisk_SK-small_scfm"),
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
  times = list(start = 0, end = 10),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_yieldTables@development",
    "PredictiveEcology/Biomass_core@development",
    "PredictiveEcology/CBM_defaults@development",
    "PredictiveEcology/CBM_dataPrep@development",
    "PredictiveEcology/CBM_core@CBM4",
    "PredictiveEcology/LandRCBM_split3pools@CBM4",
    file.path("PredictiveEcology/scfm@development/modules",
              c("scfmDataPrep",
                "scfmIgnition", "scfmEscape", "scfmSpread",
                "scfmDiagnostics")),
    "PredictiveEcology/scfm_burnProbability@main",
    "PredictiveEcology/LandRCBM_fireRisk@main"
  ),
  
  params = list(
    .globals = list(
      .plots = "png",
      .plotInterval = 10,
      sppEquivCol = "LandR",
      .studyAreaName = "SK-small",
      dataYear = 2000
    ),
    Biomass_borealDataPrep = list(
      subsetDataBiomassModel = 50
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@development",
      maxAge = 110,
      .useCache = "generateData"
    ),
    CBM_core = list(
      cohort_fixed = FALSE
    ),
    scfmDataPrep = list(
      targetN = 200,
      flammabilityThreshold = 0.05,
      .useParallelFireRegimePolys = FALSE,
      fireEpoch = c(1971, 2020),
      fireRegimePolysType = "FRU",
      # Default (6.25e8 m^2) exceeds the whole small study area's size, so every
      # fire regime polygon is flagged as a sliver with nothing left to merge
      # into; lower it well below the study area's size for this smoke test.
      sliverThreshold = 1e5
    ),
    scfm_burnProbability = list(
      # 300 iterations may be too low to generate realistic burn probability predictions.
      # nIterations = 300,
      .runInterval = 3
    ),
    LandRCBM_fireRisk = list(
      .runInterval = 1
    )
  ),
  
  packages = c("terra", "PredictiveEcology/LandR@development"),
  
  # Study area is near La Ronge, SK -- boreal forest
  studyAreaPolygon = {
    LandR::randomStudyArea(
      center = terra::vect(
        data.frame(lon = -105.3, lat = 55.1),
        geom = c("lon", "lat"),
        crs = "EPSG:4326"
      ), 
      size = 3 * 10^7, seed = 7351) |>
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
  
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea, dPath = paths$inputPath)
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv[KNN != "" & LANDIS_traits != ""] # avoid a bug with shore pine
  },
  
  disturbanceMeta = data.frame(
    eventID          = 1,
    name             = "Wildfire",
    sourceObjectName = "rstCurrentBurn",
    sourceDelay      = 1
  )
)

out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)

