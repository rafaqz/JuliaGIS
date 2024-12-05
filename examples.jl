
# Makie and Tyler.jl

using Tyler, GLMakie, Extents, NaturalEarth
using Rasters, RasterDataSources, ArchGDAL
using DataFrames
using ImageIO
import Proj
import GeometryOps as GO
import GeoInterface as GI
GLMakie.activate!()
ENV["RASTERDATASOURCES_PATH"] = ".";

countries = naturalearth("ne_10m_admin_0_countries") |> DataFrame
norway = subset(countries, :NAME => ByRow(==("Norway"))).geometry[1]
norway_webmercator = GO.resample(norway; target_crs=EPSG(3857))

# Get december mean temperatures
tmax_dec = mask(Raster(WorldClim{Climate}, :tmax; month=December, replace_missing=true); with=norway)
# Mask tmax with Norway
# Project to web mercator, which is what most web tiles use
tmax_webmercator = Rasters.reproject(Rasters.shiftlocus(Rasters.Center(), tmax_dec_masked); crs=EPSG(3857))

# Make a tiling map
tyler = Tyler.Map(GI.extent(norway); provider=Tyler.TileProviders.Google())

# Add a temperature layer
pt = Makie.plot!(tyler.axis, tmax_webmercator; alpha=0.7, colormap=:seaborn_icefire_gradient)
Makie.translate!(pt, 0, 0, 10)

# Add norway borders
p = Makie.plot!(tyler.axis, norway_webmercator; color=:transparent, strokewidth=2)
Makie.translate!(p, 0, 0, 20)

# Examples

using Rasters
using ArchGDAL
using Dates
using DataFrames
using GBIF2
using RasterDataSources
using Plots
using Rasters: trim


## `extract`: species distribution modelling

# Extract climate data at specific points:

using SpeciesDistributionModels, Rasters, RasterDataSources, GBIF2, LibGEOS, ArchGDAL, StatsBase
import GeometryOps as GO
import SpeciesDistributionModels as SDM

# We start by extracting occurrence records from GBIF using GBIF2.

# extract occurrence records for Anopheles rufipes
nili = species_match("Anopheles nili")
occurrences = occurrence_search(nili, hasCoordinate = true, year = (1990, 2020), limit = 1000)
presence_points = unique(occurrences.geometry)

# The region of interest will be defined by a convex hull around our presence points, buffered by 5 degrees.

roi = GO.buffer(GO.convex_hull(presence_points), 5)

# Now we'll load bioclimatic variables from WorldClim, using RasterDataSources and Rasters.

# load bioclimatic variables
bio = crop(RasterStack(WorldClim{BioClim}); to = roi)

# The extract function from Rasters makes it easy to extract the variables at our presence points.

# get values at presence points
presences = extract(bio, presence_points; skipmissing = true, geometry = false)

# For background (pseudo-absence) points, we draw 300 points uniformly distributed within the region of interest.

# get a raster with dimensions of bio that is true for all cells within roi and where bio is not missing
roi_raster = rasterize(roi; to = bio, fill = true, missingval = false) .* Rasters.boolmask(bio)
# draw random cells in the region of interest
background_points = sample(DimIndices(roi_raster), weights(roi_raster), 300)
# get bioclimatic variables at the background locations
background = map(p -> bio[p...], background_points)

# Let's choose some settings for modelling. Here we're using 4 different models, resample using 3-fold stratified cross-validation, and use 3 predictor variables.

# choose models, a resampler, and predictor variables
models = [SDM.linear_model(), SDM.random_forest(), SDM.random_forest(; max_depth = 3), SDM.boosted_regression_tree()]
resampler = StratifiedCV(; nfolds = 3)
predictors = (:bio1, :bio7, :bio12)

# Now we're ready to run the models. We construct an ensemble using the sdm function.

# run models and construct ensemble
ensemble = sdm(presences, background; models, resampler, predictors)

# When an ensemble is printed it shows some basic information.

SDMensemble with 12 machines across 4 groups
Occurence data: Presence-Absence with 27 presences and 300 absences 
Predictors: bio1 (Continuous), bio7 (Continuous), bio12 (Continuous)

# To see how our models performed, use the evaluate function. Here we use it with default settings.

evaluation = SDM.evaluate(ensemble)

# SDMensembleEvaluation with 4 performance measures
# train

# Finally, we predict back to a raster, taking the simple mean of all models to generate a final prediction.

pr = SDM.predict(ensemble, bio; reducer = mean)


## `zonal` statistics: find the hottest countries

using Rasters, RasterDataSources, ArchGDAL, Dates, DataFrames, NaturalEarth, Statistics

# Find the hottest and coldest countries in July:
# First get climate data
tmax = Raster(WorldClim{Climate}, :tmax; month=July)

# Then some country borders:
countries = naturalearth("ne_10m_admin_0_countries") |> DataFrame

# Now we can get some zonal statistics, and make them a new column of the DataFrame:
countries.july_maxtemp = zonal(mean, tmax; 
    of=countries, boundary=:touches, progress=false
)

# Then we can subset and sort to get the hottest:
sub = subset(countries, :july_maxtemp => ByRow(!isnan))
sort!(sub, :july_maxtemp).NAME