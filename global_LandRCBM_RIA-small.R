
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
               outputPath  = file.path(projectPath, "outputs", "LandRCBM_RIA-small"),
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
  times = list(start = 2020, end = 2022),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@development",
    "PredictiveEcology/Biomass_speciesFactorial@development",
    "PredictiveEcology/Biomass_speciesParameters@development",
    "PredictiveEcology/Biomass_yieldTables@development",
    "PredictiveEcology/Biomass_core@development",
    "PredictiveEcology/CBM_defaults@development",
    "PredictiveEcology/CBM_dataPrep@development",
    "PredictiveEcology/CBM_core@development",
    "PredictiveEcology/LandRCBM_split3pools@main"
  ),
  
  params = list(
    .globals = list(
      .plots = "png",
      .plotInterval = 10,
      sppEquivCol = "LandR",
      .studyAreaName = "RIA-small",
      dataYear = 2020
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
  
  packages = c("reproducible", "googledrive", "sf", "terra", "PredictiveEcology/LandR@development"),
  
  # Study area is a small area in the RIA
  studyAreaPolygon = {
    reproducible::prepInputs(
      url = "https://drive.google.com/file/d/1LxacDOobTrRUppamkGgVAUFIxNT4iiHU",
      destinationPath = paths$inputPath
    ) |>
      subset(TSA_NUMBER %in% c("08", "16", "24", "40", "41")) |>
      sf::st_union() |>
      sf::st_crop(c(xmin = 1000000, xmax = 1200000, ymin = 1100000, ymax = 1300000)) |>
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
  rasterToMatch_biomassParam = rasterToMatch,
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

