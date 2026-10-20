###
###
# This script runs LandR and CBM in NWT with fire simulated with scfm.
# Some data are retrieved from Google Drive, so R will ask to sign up to a Google Drive account.
# Access to permanent Sample Plot data is restricted and will require user to ask permission.
###
###

# Get the minimal amount of packages
repos <- c("predictiveecology.r-universe.dev", getOption("repos"))
if (!require("SpaDES.project")){
  Require::Install(c("SpaDES.project", "SpaDES.core", "reproducible"), repos = repos, dependencies = TRUE)
}

out <- SpaDES.project::setupProject(
  paths = list(projectPath = getwd(),
               inputPath = "inputs",
               outputPath = "outputs/SCFM",
               cachePath = "cache"),
  options = options(
    repos = c(repos = repos),
    Require.cloneFrom = Sys.getenv("R_LIBS_USER"),
    spades.moduleCodeChecks = FALSE,
    spades.recoveryMode = FALSE,
    reproducible.useMemoise = FALSE),
  times = list(start = 2020, end = 2520),
  modules = c(
    "PredictiveEcology/Biomass_borealDataPrep@0d9a34d520587156caab19de52bd184d312216b8",
    "PredictiveEcology/Biomass_speciesFactorial@dac4cc052e0e780d752faacb2250d56b68c84d11",
    "PredictiveEcology/Biomass_speciesParameters@305be17e873fcb789e681b41a2f680632ce833fc",
    "PredictiveEcology/CBM_defaults@a498ab84e243dc9b8e0d431d553fc2dae547d267",
    "PredictiveEcology/Biomass_regeneration@1d967adb38205e978329739aa110f248e70252ad",
    "PredictiveEcology/Biomass_yieldTables@a49d5ef6b1a48868b932c0030aeaa1413b4f57c1",
    "PredictiveEcology/Biomass_core@0639f3343dcf21b120df6b492422273cf81d30bd",
    "PredictiveEcology/CBM_dataPrep@f6be1956be584d9f347fd2011bf92efe1ec924f7",
    "PredictiveEcology/LandRCBM_split3pools@7359d10a33cd86cab8c013285c00112783cab888",
    "PredictiveEcology/CBM_core@ab47a77621e3088efa182b42a59bb24963015600",
    file.path("PredictiveEcology/scfm@346416fc755f50b59e29460c2e3fe1b42ab22bde/modules",
              c("scfmDataPrep",
                "scfmIgnition", "scfmEscape", "scfmSpread",
                "scfmDiagnostics"))
  ),
  outputs = data.frame(
    objectName = c("summaryBySpecies", "rasterToMatch", "disturbanceEvents"),
    file = c("summaryBySpecies.rds", "rasterToMatch.tif", "disturbanceEvents.rds"),
    fun = c("saveRDS", "writeRaster", "saveRDS"),
    package = c("base", "terra", "base")
  ),
  packages = c("googledrive", 'RCurl', 'XML', "stars", "httr2"),
  # Study area is the taiga plains of northwest territories
  studyArea = {
    targetCRS <- terra::crs(terra::rast("~/inputs/SCANFI_att_age_S_2020_v1_1.tif"))
    # northwest territories boundaries
    nwt <- reproducible::prepInputs(url = "https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/files-fichiers/lpr_000a21a_e.zip",
                                    destinationPath = "inputs",
                                    projectTo = targetCRS)
    nwt <- nwt[nwt$PRENAME == "Northwest Territories",]
    # ecozone: Taiga plains
    taigaPlains <- reproducible::prepInputs(url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
                                            destinationPath = "inputs", projectTo = nwt)
    taigaPlains <- taigaPlains[taigaPlains$ECOZONE == 4, ]
    sa <- reproducible::postProcessTo(taigaPlains, cropTo = nwt, maskTo = nwt) |>
      reproducible::Cache() |>
      sf::st_union() |>
      sf::st_as_sf()
    sa
  },
  studyArea_biomassParam = studyArea,
  rasterToMatch = {
    sa <- terra::vect(studyArea)
    rtm <- terra::rast(sa, res = c(240, 240))
    terra::crs(rtm) <- terra::crs(sa)
    sa$id <- 1
    rtm <- terra::rasterize(sa, rtm, field = "id", touches = FALSE)
    rtm
  },
  masterRaster = {
    masterRaster = rasterToMatch
    masterRaster
  },
  adminLocator = "Northwest Territories",
  ecoLocator = 4,
  sppEquiv = {
    speciesInStudy <- LandR::speciesInStudyArea(studyArea,
                                                dPath = "inputs")
    species <- LandR::equivalentName(speciesInStudy$speciesList, df = LandR::sppEquivalencies_CA, "LandR")
    sppEquiv <- LandR::sppEquivalencies_CA[LandR %in% species]
    sppEquiv <- sppEquiv[KNN != "" & LANDIS_traits != ""] #avoid a bug with shore pine
  },
  disturbanceMeta = data.table(eventID = 1,
                               disturbance_type_id = 1,
                               wholeStand = 1,
                               name = "Wildfire",
                               sourceValue = 1,
                               sourceDelay = 1,
                               sourceObjectName = "rstCurrentBurn"),
  params = list(
    .globals = list(
      .plots = c("png"),
      .plotInterval = 10,
      sppEquivCol = 'LandR',
      .studyAreaName = "NWT",
      dataYear = 2020
    ),
    CBM_core = list(
      skipPrepareCBMvars = TRUE,
      .saveInterval = 10
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
      .plots = "png",
      standAgesForFitting = c(0, 100),
      .useCache = c(".inputObjects", "init"),
      speciesFittingApproach = "focal"
    ),
    Biomass_yieldTables = list(
      moduleNameAndBranch = "PredictiveEcology/Biomass_core@0639f3343dcf21b120df6b492422273cf81d30bd",
      maxAge = 150,
      .plots = "png",
      .useCache = "generateData"
    ),
    scfmDataPrep = list(targetN = 4000,
                        flammabilityThreshold = 0.05,
                        .useParallelFireRegimePolys = FALSE,
                        fireEpoch = c(1971, 2020),
                        fireRegimePolysType = "FRU"
    )
  )
)

out$loadOrder <- unlist(out$modules)

initOut <- SpaDES.core::simInit2(out)
simOut <- SpaDES.core::spades(initOut)

