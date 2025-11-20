
## TODO: SDA aggregation may contain some values that are too low
# https://github.com/ncss-tech/soilDB/issues/257

## get the latest
# remotes::install_github("brownag/rgeedim")

## comparison ideas
# https://gist.github.com/brownag/edabb80e8b9cbc2cb3057e998556e08c

## modified SOLUS app
# https://code.earthengine.google.com/e264b766650362730527e09756c0ef2e
library(reticulate)
use_python("C:/Python310/python.exe", required=T)
ee = import("ee")
ee$Initialize()


library(rgeedim)
library(soilDB)
library(aqp)
library(terra)
library(viridis)
library(rasterVis)
library(latticeExtra)
library(tactile)


source('./code/SOLUS/local-functions.R')

txt.lut <- read.csv('http://soilmap2-1.lawr.ucdavis.edu/800m_grids/RAT/texture_2550.csv')

# this should bring up a web page with steps to setup temporary auth
# use the default EE project
# TODO: document this!
ee = import("ee")
# use as needed
gd_authenticate(auth_mode="notebook")

# init python back-end
gd_initialize()

## TODO: create a list of variables in each data source and possibly required rescaling factor

# variable of interest
.variable <- 'awc_r'
.SOLUS_variable <- 'AWC'
.rescale_value <- 0.01

.variable <- 'claytotal_r'
.SOLUS_variable <- 'Clay'
.rescale_value <- 1

.variable <- 'sandtotal_r'
.SOLUS_variable <- 'Sand'
.rescale_value <- 1

.variable <- 'ec_r'
.SOLUS_variable <- 'EC'
.rescale_value <- 0.1

.variable <- 'ph1to1h2o_r'
.SOLUS_variable <- 'pH'
.rescale_value <- 0.01




## SoilWeb style bounding boxes

# MN examples from Joe Brennan
bb <- '-96.4939 46.4092,-96.4939 46.6174,-95.8348 46.6174,-95.8348 46.4092,-96.4939 46.4092'



# TX027
bb <- '-97.4535 30.8660,-97.4535 31.1488,-96.9279 31.1488,-96.9279 30.8660,-97.4535 30.8660'
# zoom
bb <- '-96.7651 30.8599,-96.7651 30.9307,-96.6337 30.9307,-96.6337 30.8599,-96.7651 30.8599'


# TX331 | TX395
bb <- '-97.0010 30.7474,-97.0010 31.0306,-96.4754 31.0306,-96.4754 30.7474,-97.0010 30.7474'

# Capitol Reef
# there should be abrupt changes in AWC
bb <- '-111.3555 38.2650,-111.3555 38.3944,-111.0927 38.3944,-111.0927 38.2650,-111.3555 38.2650'


# TX
bb <- '-97.3052 30.5714,-97.3052 30.8551,-96.7796 30.8551,-96.7796 30.5714,-97.3052 30.5714'


# CA630
bb <- '-121 37.9901,-121 38.1200,-120.2 38.1200,-120.2 37.9901,-121 37.9901'

# sand hills, southern interface
bb <- '-101.4416 41,-101.4416 41.3678,-100.9159 41.3678,-100.9159 41,-101.4416 41'


# Tulare lake basin
bb <- '-120.2591 35.9366,-120.2591 36.2033,-119.7335 36.2033,-119.7335 35.9366,-120.2591 35.9366'


# western Fresno county
bb <- '-120.6989 36.4359,-120.6989 36.7009,-120.1733 36.7009,-120.1733 36.4359,-120.6989 36.4359'


## check assumptions on misc. areas etc.
# high sierra
bb <- '-119.6 36.1152,-119.6 36.8857,-118.5 36.8857,-118.5 36.1152,-119.6 36.1152'


# RI600
bb <- '-71.8018 41.3583,-71.8018 41.6054,-71.2762 41.6054,-71.2762 41.3583,-71.8018 41.3583'

# Valley Springs
bb <- '-121.0269 38.0948,-121.0269 38.2246,-120.7641 38.2246,-120.7641 38.0948,-121.0269 38.0948'

# near Chesapeak Bay
bb <- '-76.0056 38.7209,-76.0056 38.8495,-75.7428 38.8495,-75.7428 38.7209,-76.0056 38.7209'


# Mississippi river
bb <- '-90.4532 33.5122,-90.4532 33.6496,-90.1904 33.6496,-90.1904 33.5122,-90.4532 33.5122'
bb <- '-90.5847 33.4435,-90.5847 33.7183,-90.0591 33.7183,-90.0591 33.4435,-90.5847 33.4435'

# from Chad, coastal plain
bb <- '-80.3176 33.9643,-80.3176 34.0977,-79.9880 34.0977,-79.9880 33.9643,-80.3176 33.9643'


# favorite landscapes
bb <- '-90.1898 41.8027,-90.1898 41.8642,-90.0584 41.8642,-90.0584 41.8027,-90.1898 41.8027'

# Red Barn
bb <- '-119.7765 36.8951,-119.7765 36.9611,-119.6451 36.9611,-119.6451 36.8951,-119.7765 36.8951'


# LA, Gulf interface
# bb <- '-91.8066 29.4575,-91.8066 30.0305,-90.7553 30.0305,-90.7553 29.4575,-91.8066 29.4575'

# Pinnacles National Park
bb <- '-121.3049 36.4195,-121.3049 36.5521,-121.0420 36.5521,-121.0420 36.4195,-121.3049 36.4195'

# Salinas Valley
bb <- '-121.6715 36.4873,-121.6715 36.6198,-121.4087 36.6198,-121.4087 36.4873,-121.6715 36.4873'

# Kings River alluvial fan / outwash sequences
bb <- '-119.7997 36.6051,-119.7997 36.8695,-119.2741 36.8695,-119.2741 36.6051,-119.7997 36.6051'

# Gabbro + vertisols near Sanger, CA
bb <- '-119.5323 36.6515,-119.5323 36.7837,-119.2695 36.7837,-119.2695 36.6515,-119.5323 36.6515'


# Turlock Lake, CA
bb <- '-120.7047 37.5502,-120.7047 37.6808,-120.4419 37.6808,-120.4419 37.5502,-120.7047 37.5502'

# Chico, CA
bb <- '-121.9556 39.6609,-121.9556 39.7878,-121.6928 39.7878,-121.6928 39.6609,-121.9556 39.6609'

# Sacramento River, Glenn | Butte co. boundary
bb <- '-122.1354 39.5018,-122.1354 39.6290,-121.8725 39.6290,-121.8725 39.5018,-122.1354 39.5018'

# Putah Creek
bb <- '-122.2131 38.6740,-122.2131 38.8027,-121.9503 38.8027,-121.9503 38.6740,-122.2131 38.6740'




# create WKT from SoilWeb style BBOX
bb <- sprintf("POLYGON((%s))", bb)



tx.clay <- getData(bb, .variable = 'claytotal_r', .SOLUS_variable = 'Clay', .rescale_value = 1)
tx.sand <- getData(bb, .variable = 'sandtotal_r', .SOLUS_variable = 'Sand', .rescale_value = 1)

g.texture <- tx.sand$gNATSGO
values(g.texture) <- ssc_to_texcl(sand = values(tx.sand$gNATSGO), clay = values(tx.clay$gNATSGO), droplevels = TRUE)

s.texture <- tx.sand$gNATSGO
values(s.texture) <- ssc_to_texcl(sand = values(tx.sand$SOLUS), clay = values(tx.clay$SOLUS), droplevels = TRUE)

## match colors
g.cols <- txt.lut$hex[match(
  levels(g.texture)[[1]]$label,
  txt.lut$class
)]

s.cols <- txt.lut$hex[match(
  levels(s.texture)[[1]]$label,
  txt.lut$class
)]

ragg::agg_png(filename = 'PINN-texture-15cm.png', width = 2800, height = 1250, scaling = 2)

par(mar = c(0, 0, 0, 0), mfrow = c(1, 2))
plot(g.texture, axes = FALSE, col = g.cols, main = 'gNATSGO\nTexture Class <2mm\n15cm', line = -2)
plot(s.texture, axes = FALSE, col = s.cols, main = 'SOLUS\nTexture Class <2mm\n15cm', line = -2)

dev.off()

# check
levelplot(c(tx.clay$gNATSGO, tx.clay$SOLUS),
          scales = list(alternating = 1),
          col.regions = magma,
          main = sprintf('%s near 15cm', 'Clay'),
          names.attr = c('gNATSGO', 'SOLUS'),
          sub = 'EPSG 5070, 30m grid cell size'
)



