# type your code here
library(raster)



clim_list <- list.files("data/bioclim/",pattern=".bil$",full.names = T) # '..' leads to the path above the folder where the .rmd file is located

bio1 <- clim[[1]]

length(values(bio1))
nrow(bio1) * ncol(bio1)

values(bio1)  [values(bio1) > 0]  <- -100
plot(bio1)


View(occ_raw)


occ_raw[1,]
names(occ_raw)
plot(occ_unique)
myCRS1 <- CRS("+init=epsg:4326")

crs(occ_unique)
