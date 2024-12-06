nothing

# Makie and Tyler.jl example
using Tyler # Tiled plotting in Makie.jl
using GLMakie # Open gl plotting backend
using Extents # Shared Extent object
using NaturalEarth # Courntry borders
using Rasters # Raster data
using RasterDataSources # Climate data
using ArchGDAL # GDAL bindings
using DataFrames # DataFrames
using Dates # Date and time handling
using StatsBase
import Proj # Projections
import GeometryOps as GO # Spatial operations
import GeoInterface as GI # Common interface for geospatial objects
ENV["RASTERDATASOURCES_PATH"] = "."; # Just for CI

# Get Norway country borders
countries = naturalearth("ne_10m_admin_0_countries") |> DataFrame
norway = subset(countries, :NAME => ByRow(==("Norway"))).geometry[1]
norway_webmercator = GO.reproject(norway; target_crs=EPSG(3857))

# Get december mean temperatures 
tmax_dec = Raster(WorldClim{Climate}, :tmax; month=December, replace_missing=true)

# Mask with Norway
tmax_dec_masked = mask(tmax_dec; with=norway, boundary=:touches)

# Project to web mercator, which is what most web tiles use
tmax_webmercator = Rasters.reproject(Rasters.shiftlocus(Rasters.Center(), tmax_dec_masked); crs=EPSG(3857))

# Make a tiling map
tyler = Tyler.Map(GI.extent(norway); provider=Tyler.TileProviders.Google())

# Add a temperature layer
pt = Makie.plot!(tyler.axis, tmax_webmercator; alpha=0.9, colormap=:seaborn_icefire_gradient)
Makie.translate!(pt, 0, 0, 10)

# Add norway borders
p = Makie.plot!(tyler.axis, norway_webmercator; color=:transparent, strokewidth=2)
Makie.translate!(p, 0, 0, 20)


# SDM Example

using Rasters, RasterDataSources, ArchGDAL, NaturalEarth, DataFrames
bio = RasterStack(WorldClim{BioClim}, (1,12))
countries = naturalearth("ne_10m_admin_0_countries") |> DataFrame

# Subset australia from the data
australia = subset(countries, :NAME => ByRow(==("Australia"))).geometry
bio_aus = Rasters.trim(mask(bio; with = australia)[X = 110 .. 156, Y = -45 .. -10])

# Let's plot this data to see what it looks like.
using CairoMakie
CairoMakie.activate!()
Rasters.rplot(bio_aus)

## Occurrence data
# Next, we use [GBIF2.jl](www.github.com/rafaqz/GBIF2.jl) to download occurrence records for this species. We use the [thin](@ref) function in this package to weed out occurrences that are very close to each other, using a cut-off of 5km.
using GBIF2, SpeciesDistributionModels
sp = species_match("Eucalyptus regnans")
occurrences_raw = occurrence_search(sp; year = (1970,2000), country = "AU", hasCoordinate = true, limit = 2000)
occurrences = thin(occurrences_raw.geometry, 5000)

## Background points
# Next, we sample random points to use as background points.

# Let's plot both the occurrence and background points to see where _Eucalyptus regnans_ is found.

using StatsBase
bg_indices = sample(findall(boolmask(bio_aus)), 500)
bg_points = DimPoints(bio_aus)[bg_indices]
fig, ax, pl = Makie.plot(bio_aus.bio1)
Makie.scatter!(ax, occurrences; color = :red)
Makie.scatter!(ax, bg_points; color = :grey)
fig

## Handling data
# SpeciesDistributionModels.jl has a [sdmdata](@ref) function to handle input data. It takes tabular presence and background data as inputs, such as what is returned by `Rasters.extract` and `Rasters.sample`.

using SpeciesDistributionModels
p_data = extract(bio_aus, occurrences; skipmissing = true)
bg_data = bio_aus[bg_indices]
data = sdmdata(p_data, bg_data; resampler = CV(nfolds = 3))

## Fitting an ensemble
# Now that we have our `data` object with presence and background data, we can fit our ensemble. The `sdm` function fits a whole ensemble, taking two arguments: a data object and a `NamedTuple` with models the ensemble should have. This can be any MLJ-compatible model. In this case, we use Maxnet, boosted regression trees (from the EvoTrees.jl package), and a GLM.

using Maxnet: MaxnetBinaryClassifier
using EvoTrees: EvoTreeClassifier
using MLJGLMInterface: LinearBinaryClassifier
models = (
    maxnet = MaxnetBinaryClassifier(),
    brt = EvoTreeClassifier(),
    glm = LinearBinaryClassifier()
)

ensemble = sdm(data, models)

## Evaluating an ensemble
import SpeciesDistributionModels as SDM
ev = SDM.evaluate(ensemble; measures = (; auc, accuracy))


## Predicting
pred = SDM.predict(ensemble, bio_aus; reducer = mean)
Makie.plot(pred; colorrange = (0, 1))


## Understanding the model
expl = SDM.explain(ensemble; method = ShapleyValues(8))
variable_importance(expl)


# Interactively plot the model explanation
using GLMakie
GLMakie.activate!()
interactive_response_curves(expl)





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