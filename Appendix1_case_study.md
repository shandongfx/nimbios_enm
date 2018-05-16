---
title: "A brief tutorial on running Maxent in R (v1.1)"
author: "Xiao Feng, Cassondra Walker, and Fikirte Gebresenbet"
date: "May 16, 2018"
output:
  html_document:
    keep_md: yes
    toc: yes
    toc_depth: 4
  pdf_document:
    toc: yes
    toc_depth: '4'
---
##1. Set up the working environment  
###1.1 Load packages  
Running Maxent in R requires several packages. Specifically, the "dismo" package, which contains *maxent* function that calls *maxent.jar* in R, the *raster* package, which provides functions for analyzing gridded data, the *rgeos* package, which provides functions for analyzing spatial data.

#####Thread 1

```r
packages_needed <- c("raster", # for raster analysis
                     "dismo", # a collection of ENM/SDM tools
                     "rgeos","rgdal","sp", # spatial data analysis
                     "ENMeval", # a few new tools in ENM/SDM
                     "utils", # for zip & unzip files
                     "jsonlite" # necessary for download data from GBIF
                     )
pk_to_install <- packages_needed [!( packages_needed %in% rownames(installed.packages())  )]
if(length(pk_to_install)>0 ){
  install.packages(pk_to_install,repos="http://cran.r-project.org")
}
#lapply(packages_needed, require, character.only = TRUE)
library("raster")
library("dismo")
library("rgeos")
library("rgdal")
library("sp")
library("ENMeval")
```
																												   
#####Thread 2

```r
if( !("rJava" %in% rownames(installed.packages()))  ){
  install.packages("rJava",repos="http://cran.r-project.org")
}
#If you are using a Mac machine, an additional step may be needed before loading rJava package
#dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_144.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_171.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
library("rJava")
```




###1.2 Set up the Maxent path  
In order for Maxent to work properly in R, the *maxent.jar* file needs to be accessible by *dismo* package.  

#####Thread 4

```r
# download maxent.jar 3.3.3k, and place the file in the desired folder; note that, there may be a newer version of Maxent
if( !file.exists(paste0(system.file("java", package="dismo"),"/maxent.jar"))  )   {
utils::download.file(url="https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar",
                     destfile=paste0(system.file("java", package="dismo"),"/maxent.jar"),
                     mode="wb") ## wb for binary file, otherwise maxent.jar can not execute
}
# also note that both R and Java need to be the same bit (either 32 or 64) to be compatible to run

# to increase memory size of the JVM and prevent memory issues with Maxent.jar
# options( java.parameters = c("-Xss2560k", "-Xmx2g") ) 
```


##2. Prepare data input  
###2.1 Load environmental layers 
In our example, we use bioclimatic variables (downloaded from worldclim.org) as input environmental layers. We stack our environmental layers so that they can be processed simultaneously.  

#####Thread 5 

```r
# prepare folders for data input and output
if(!file.exists("data")) dir.create("data")
if(!file.exists("data/bioclim")) dir.create("data/bioclim")
if(!file.exists("data/studyarea")) dir.create("data/studyarea")
if(!file.exists("output")) dir.create("output")
require(utils)
# download climate data from worldclim.org
if( !file.exists( paste0("data/bioclim/bio_10m_bil.zip")   )){
utils::download.file(url="http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip",
                     destfile=paste0("data/bioclim/bio_10m_bil.zip")) 
utils::unzip("data/bioclim/bio_10m_bil.zip",exdir="data/bioclim/")
}
# This searches for all files that are in the path "data/bioclim/" and 
# have a file extension of .bil. You can edit this code to reflect the path name 
# and file extension for your environmental variables
clim_list <- list.files("data/bioclim/",pattern=".bil$",full.names = T) # '..' leads to the path above the folder where the .rmd file is located

# stacking the bioclim variables to process them at one go 
clim <- raster::stack(clim_list) 
```

####2.1.1 simple Raster manipulation    

```r
# we want the new layer to be 10 times coarser at each axis (100 times coarser)
# read current resolution
bio1 <- clim[[1]]
bio1
```

```
## class       : RasterLayer 
## dimensions  : 900, 2160, 1944000  (nrow, ncol, ncell)
## resolution  : 0.1666667, 0.1666667  (x, y)
## extent      : -180, 180, -60, 90  (xmin, xmax, ymin, ymax)
## coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
## data source : /Users/labmac/Desktop/ENM_demo/peerj example/workshop_maxent_R-master/data/bioclim/bio1.bil 
## names       : bio1 
## values      : -269, 314  (min, max)
```

```r
nrow(bio1)
```

```
## [1] 900
```

```r
ncol(bio1)
```

```
## [1] 2160
```

```r
extent(bio1)
```

```
## class       : Extent 
## xmin        : -180 
## xmax        : 180 
## ymin        : -60 
## ymax        : 90
```

```r
# define new resolution
newRaster <- raster( nrow= nrow(bio1)/10 , ncol= ncol(bio1)/10 )

# define extent
extent(newRaster) <- extent(bio1)

# fill the new layer with new values
newRaster <- resample(x=bio1,y=newRaster,method='bilinear')
plot(bio1)
```

![](Appendix1_case_study_files/figure-html/raster_manipulation-1.png)<!-- -->

```r
plot(newRaster) # new layer seems coarser
```

![](Appendix1_case_study_files/figure-html/raster_manipulation-2.png)<!-- -->

####2.1.2 reclassify raster layer    


```r
# we want to classify the world into two classes based on temperature, high > mean & low < mean
myLayer<- raster("data/bioclim/bio1.bil")

# values smaller than meanT becomes 1; values larger than meanT will be 2
myMethod <- c(-Inf, 100, 1,  100, Inf, 2)
myLayer_classified <- reclassify(myLayer,rcl= myMethod)
plot(myLayer_classified)
```

![](Appendix1_case_study_files/figure-html/reclassify_raster-1.png)<!-- -->

####2.1.3 raster calculation    


```r
# read precipitation data 
wet <- raster("data/bioclim/bio13.bil") # precipitation of wettest month
dry <- raster("data/bioclim/bio14.bil") # precipitation of driest month
plot(stack(wet,dry))
```

![](Appendix1_case_study_files/figure-html/raster_calculation-1.png)<!-- -->

```r
# calculate the difference
diff <- wet - dry
#plot(diff)

# calculate the mean of the two month
twoLayers <- stack(wet,dry)
meanPPT <- calc(twoLayers,fun=mean)
#plot(meanPPT)

# the following code gives the same results
meanPPT2 <-  (wet + dry)/2

# histogram of one layer
hist(twoLayers[[1]])
```

```
## Warning in .hist1(x, maxpixels = maxpixels, main = main, plot = plot, ...):
## 5% of the raster cells were used. 100000 values used.
```

![](Appendix1_case_study_files/figure-html/raster_calculation-2.png)<!-- -->

```r
# correlation between different layers
pairs(twoLayers[[1:2]])
```

![](Appendix1_case_study_files/figure-html/raster_calculation-3.png)<!-- -->


###2.2 Occurrence data  
####2.2.1 Download occurrence data  
For our example, we download occurrence data of the nine-banded armadillo from GBIF.org (Global Biodiversity Information Facility). 

#####Thread 6

```r
# We provide an if/else statement that checks for occurrence data that have already been downloaded.
# The goal of this module is to avoid downloading the occurrence data multiple times.
require(jsonlite)				 
```

```
## Loading required package: jsonlite
```

```r
if(file.exists("data/occ_raw")){
  load("data/occ_raw")
}else{
  occ_raw <- gbif(genus="Dasypus",species="novemcinctus")
  save(occ_raw,file = "data/occ_raw")
  write.csv("data/occ_raw.csv")
}

# to view the first few lines of the occurrence dataset
# head( occ_raw )
```

####2.2.2 Clean occurrence data
Since some of our records do not have appropriate coordinates and some have missing locational data, we need to remove them from our dataset. To do this, we creatd a new dataset named “occ_clean”, which is a subset of the “occ_raw” dataset where records with missing latitude and/or longitude are removed. This particular piece of code also returns the number of records that are removed from the dataset. Additionally, we remove duplicate records and create a subset of the cleaned data with the duplicates removed. 

#####Thread 7: remove data without coordinate, remove duplicated coordinates, only keep SPECIMEN data, limit the temporal range of data

```r
# remove erroneous coordinates, where either the latitude or longitude is missing
occ_clean <- subset(occ_raw,(!is.na(lat))&(!is.na(lon))) #  "!" means the opposite logic value
cat(nrow(occ_raw)-nrow(occ_clean), "records are removed")
```

```
## 2337 records are removed
```

```r
# remove duplicated data based on latitude and longitude
dups <- duplicated(occ_clean[c("lat","lon")])
occ_unique <- occ_clean[!dups,]
cat(nrow(occ_clean)-nrow(occ_unique), "records are removed")
```

```
## 1205 records are removed
```

```r
# only keep specimen
table(occ_unique$basisOfRecord)
```

```
## 
##     FOSSIL_SPECIMEN   HUMAN_OBSERVATION          LITERATURE 
##                  13                1908                   8 
## MACHINE_OBSERVATION         OBSERVATION  PRESERVED_SPECIMEN 
##                  26                  34                 906 
##             UNKNOWN 
##                 247
```

```r
occ_unique <- subset(occ_unique, basisOfRecord=="PRESERVED_SPECIMEN")
table(occ_unique$basisOfRecord)
```

```
## 
## PRESERVED_SPECIMEN 
##                906
```

```r
# limit species by year
#table(occ_unique$year)
hist(occ_unique$year)
```

![](Appendix1_case_study_files/figure-html/clean_data1-1.png)<!-- -->

```r
occ_unique <- subset(occ_unique, year>=1950)
```

Up to this point we have been working with a data frame, but it has no spatial relationship with environmental layers. So we need to make the data spatial. Once our data is spatial we can use the *plot* function to see the occurrence data and allow us to check for data points that appear to be erroneous.

#####Thread 8: make occ spatial, assign coordinate reference system to *spatial points*

```r
# make occ spatial
coordinates(occ_unique) <- ~ lon + lat

# add CRS projection.

myCRS1 <- CRS("+init=epsg:4326") # WGS 84
#myCRS2 <- CRS("+init=epsg:4269") # NAD 83
#myCRS3 <- CRS("+init=epsg:3857") # Mercator
crs(occ_unique) <- myCRS1
plot(occ_unique)
```

![](Appendix1_case_study_files/figure-html/clean_data2-1.png)<!-- -->

```r
# full reference list can be found here http://spatialreference.org/ref/

## look for erroneous points
plot(clim[[1]]) # to the first layer of the bioclim layers as a reference
plot(occ_unique,add=TRUE) # plot the oc_unique on the above raster layer
```

![](Appendix1_case_study_files/figure-html/clean_data2-2.png)<!-- -->

In the figure above, we can see several points that appear outside the known distribution of _Dasypus novemcinctus_ (North and South America) and we need to remove these from our occurrence dataset. To do this, we only kept points that have longitudes between -110° and -40°. 

#####Thread 9: remove spatial outliers, extract environmental conditions, select points by area

```r
# remove erroneous points (i.e., only keep good records)
occ_unique <- occ_unique[which(occ_unique$lon>-110 &
                                 occ_unique$lon < -40),]

# remove points without environmental data (e.g., points fall in the ocean)
conditions_occ <- extract(clim,occ_unique)
head(conditions_occ)
```

```
##      bio1 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 bio2
## [1,]  190   277    93  1049   116    59    18   306   217   217   239  131
## [2,]  156   158   154  2061   274    85    32   678   321   552   678   98
## [3,]  212   215   207  2528   353    88    37   855   347   369   855   99
## [4,]  147   148   144  1787   246    78    34   587   294   487   587   91
## [5,]  147   148   144  1787   246    78    34   587   294   487   587   91
## [6,]  190   277    93  1049   116    59    18   306   217   217   239  131
##      bio3 bio4 bio5 bio6 bio7 bio8 bio9
## [1,]   38 7101  354   18  336  229  277
## [2,]   92  172  210  104  106  154  155
## [3,]   89  322  272  161  111  207  214
## [4,]   94  180  196  100   96  144  147
## [5,]   94  180  196  100   96  144  147
## [6,]   38 7101  354   18  336  229  277
```

```r
bad_records <- is.na( conditions_occ[,1] ) 
table(bad_records)
```

```
## bad_records
## FALSE  TRUE 
##   656     2
```

```r
conditions_occ[bad_records,]
```

```
##      bio1 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 bio2
## [1,]   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA   NA
## [2,]   NA    NA    NA    NA    NA    NA    NA    NA    NA    NA    NA   NA
##      bio3 bio4 bio5 bio6 bio7 bio8 bio9
## [1,]   NA   NA   NA   NA   NA   NA   NA
## [2,]   NA   NA   NA   NA   NA   NA   NA
```

```r
occ_unique <- occ_unique[!bad_records,]

# only select points in US
country_shp <- shapefile("data/GIS polygon/country.shp")
head(country_shp)
```

```
##   FIPS_CNTRY GMI_CNTRY          CNTRY_NAME           SOVEREIGN POP_CNTRY
## 0         AA       ABW               Aruba         Netherlands     67074
## 1         AC       ATG Antigua and Barbuda Antigua and Barbuda     65212
## 2         AF       AFG         Afghanistan         Afghanistan  17250390
## 3         AG       DZA             Algeria             Algeria  27459230
## 4         AJ       AZE          Azerbaijan          Azerbaijan   5487866
## 5         AL       ALB             Albania             Albania   3416945
##    SQKM_CNTRY SQMI_CNTRY CURR_TYPE CURR_CODE LANDLOCKED COLOR_MAP
## 0     182.926     70.628    Florin       AWG          N         1
## 1     462.378    178.524 EC Dollar       XCD          N         2
## 2  641869.188 247825.703   Afghani       AFA          Y         3
## 3 2320972.000 896127.312     Dinar       DZD          N         3
## 4   85808.203  33130.551     Manat      <NA>          Y         4
## 5   28754.500  11102.110       Lek       ALL          N         6
```

```r
US_shp <- subset(country_shp,CNTRY_NAME=="United States") 
plot(US_shp)
```

![](Appendix1_case_study_files/figure-html/clean_data3-1.png)<!-- -->

```r
# this will lead to error, because of different CRS 
# over(x = occ_unique,y = US_shp) 
crs(US_shp)
```

```
## CRS arguments:
##  +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0
```

```r
crs(occ_unique)
```

```
## CRS arguments:
##  +init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
## +towgs84=0,0,0
```

```r
US_shp_prj <- spTransform(US_shp,crs(occ_unique))
over_results <- over(occ_unique,US_shp_prj)
in_us <- subset(occ_unique, !is.na(over_results$FIPS_CNTRY))
out_us <- subset(occ_unique, is.na(over_results$FIPS_CNTRY))
plot(in_us)
```

![](Appendix1_case_study_files/figure-html/clean_data3-2.png)<!-- -->

```r
plot(out_us)
```

![](Appendix1_case_study_files/figure-html/clean_data3-3.png)<!-- -->

We want to use only one occurrence point per pixel, so we need to thin our occurrence data.

#####Thread 10: spatial filter (1 point per window)

```r
# thin occ data (keep one occurrence point per cell)
cells <- cellFromXY(clim[[1]],occ_unique)
dups <- duplicated(cells)
occ_final <- occ_unique[!dups,]
cat(nrow(occ_unique)-nrow(occ_final), "records are removed")
```

```
## 122 records are removed
```

```r
# plot the first climatic layer (or replace [[1]] with any nth number of the layer of interest from the raster stack).

plot(clim[[1]]) 

# plot the final occurrence data on the environmental layer
plot(occ_final,add=T,col="red") # the 'add=T' tells R to put the incoming data on the existing layer
```

![](Appendix1_case_study_files/figure-html/clean_data4-1.png)<!-- -->

```r
# export as shapefile
shapefile(occ_final,filename="data/occ_final.shp",overwrite=TRUE)
```

```
## Warning in rgdal::writeOGR(x, filename, layer, driver = "ESRI Shapefile", :
## Field names abbreviated for ESRI Shapefile driver
```

```
## Warning in rgdal::writeOGR(x, filename, layer, driver = "ESRI Shapefile", :
## field names ht://_/, h://_/D changed by driver to: ht_//_/, h_//_/D
```

```r
# load a shapefile
# occ_final <- shapefile(filename="data/occ_final.shp")
```

###2.3 Set up study area
We create a buffer around our occurrence locations and define this as our study region, which will allow us to avoid sampling from a broad background. We establish a four-decimal-degree buffer around the occurrence points. To make sure that our buffer encompasses the appropriate area, we plot the occurrence points, the first environmental layer, and the buffer polygon.

#####Thread 11: build a buffer of occ

```r
# this creates a 4-decimal-degree buffer around the occurrence data 
occ_buff <- buffer(occ_final,width=400000) #unit is meter 

# plot the first element ([[1]]) in the raster stack
plot(clim[[1]]) 

plot(occ_final,add=T,col="red") # adds occurrence data to the plot
plot(occ_buff,add=T,col="blue") # adds buffer polygon to the plot
```

![](Appendix1_case_study_files/figure-html/set_up_study_area1-1.png)<!-- -->

With a defined study area and the environmental layers stacked, we then clip the layers to the extent of our study area. However, for ease of processing, we do this in two steps rather than one. First, we create a coarse rectangular shaped study area around the study area to reduce the size of environmental data and then extract by *mask* using the buffer we create to more accurately clip environmental layers. This approach could be faster than directly masking. We save the cropped environmental layers as .asc (ascii files) as inputs for Maxent.

#####Thread 12: crop raster by polygon (shapefile); export a group of raster layers

```r
# crop study area to a manageable extent (rectangle shaped)
studyArea <- crop(clim,extent(occ_buff))  

# the 'study area' created by extracting the buffer area from the raster stack
studyArea <- mask(studyArea,occ_buff)
# output will still be a raster stack, just of the study area

# save the new study area rasters as ascii
writeRaster(studyArea,
            # a series of names for output files
            filename=paste0("data/studyarea/",names(studyArea),".asc"), 
            format="ascii", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)
```

In the next step, we selected 10,000 random background points from the study area. To make our experiment reproducible (i.e., select the same set of points), we used a static seed via *set.seed(1)* function. Then, we plotted the background points together with the study area and occurrence data.

#####Thread 13: manually select random bg points

```r
# select background points from this buffered area; when the number provided 
# to set.seed() function, the same random sample will be selected in the next line			
# use this code before the sampleRandom function every time, if you want to get
# the same "random samples"
set.seed(1) 
bg <- sampleRandom(x=studyArea,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 

plot(studyArea[[1]])
# add the background points to the plotted raster
plot(bg,add=T) 
# add the occurrence data to the plotted raster
plot(occ_final,add=T,col="red")
```

![](Appendix1_case_study_files/figure-html/set_up_study_area3-1.png)<!-- -->

###2.4 Split occurrence data into training & testing
We randomly selected 50% of the occurrence data for model training and used the remaining for model testing. To make our experiment reproducible (i.e., select the same set of points), we used a static seed via *set.seed(1)* function.

#####Thread 14: random cut, block/spatial cut

```r
# get the same random sample for training and testing
set.seed(1) 

# randomly select 50% for training
selected <- sample(1:nrow(occ_final),nrow(occ_final)*0.5)

occ_train <- occ_final[selected,] # this is the selection to be used for model training
occ_test <- occ_final[-selected,] # this is the opposite of the selection which will be used for model testing
plot(occ_train,col="blue")
plot(occ_test,col="red",add=T)
```

![](Appendix1_case_study_files/figure-html/cut_occ_into_training_testing-1.png)<!-- -->

```r
### add another method, block cut
#library(ENMeval)
cut_block <- ENMeval::get.block(occ=as.data.frame(occ_final@coords), 
                       bg.coords=as.data.frame(bg@coords))
occ_final@data$cut_block <- cut_block$occ.grp
bg@data$cut_block <- cut_block$bg.grp

head(occ_final)
```

```
##      acceptedNameUsage
## 1413              <NA>
## 1614              <NA>
## 1701              <NA>
## 1716              <NA>
## 1819              <NA>
## 1958              <NA>
##                                                                                                                           accessRights
## 1413 Open Access, http://creativecommons.org/publicdomain/zero/1.0/; see Yale Peabody policies at: http://hdl.handle.net/10079/8931zqj
## 1614                                                                                                                              <NA>
## 1701                                                                                                                              <NA>
## 1716                                                                                                                              <NA>
## 1819                                                                                           http://vertnet.org/resources/norms.html
## 1958                                                                                           http://vertnet.org/resources/norms.html
##          adm1         adm2 associatedReferences      basisOfRecord
## 1413    Texas  Leon County                 <NA> PRESERVED_SPECIMEN
## 1614   Tolima     Planadas                 <NA> PRESERVED_SPECIMEN
## 1701   Tolima    Chaparral                 <NA> PRESERVED_SPECIMEN
## 1716   Tolima    Chaparral                 <NA> PRESERVED_SPECIMEN
## 1819 Oklahoma Noble County                 <NA> PRESERVED_SPECIMEN
## 1958 Missouri      Hickory                 <NA> PRESERVED_SPECIMEN
##      behavior                           bibliographicCitation
## 1413     <NA> Dasypus novemcinctus mexicanus (YPM MAM 015982)
## 1614     <NA>                                            <NA>
## 1701     <NA>                                            <NA>
## 1716     <NA>                                            <NA>
## 1819     <NA>                                            <NA>
## 1958     <NA>                                            <NA>
##       catalogNumber    class classKey
## 1413 YPM MAM 015982 Mammalia      359
## 1614     CZUT-M1764 Mammalia      359
## 1701     CZUT-M1766 Mammalia      359
## 1716     CZUT-M1770 Mammalia      359
## 1819          37172 Mammalia      359
## 1958         693990 Mammalia      359
##                                                                                          cloc
## 1413                             Dunn Ranch, Texas, Leon County, United States, NORTH_AMERICA
## 1614                             Vereda San Miguel, Tolima, Planadas, Colombia, SOUTH_AMERICA
## 1701                         Vereda Vega Chiquita, Tolima, Chaparral, Colombia, SOUTH_AMERICA
## 1716                           Vereda La Virginia, Tolima, Chaparral, Colombia, SOUTH_AMERICA
## 1819 0.6 mi S Perry exit, Interstate 35, Oklahoma, Noble County, United States, NORTH_AMERICA
## 1958       Route 6, near Jones Spring, Avery, Missouri, Hickory, United States, NORTH_AMERICA
##        collectionCode                    collectionID     continent
## 1413               VZ                            <NA> NORTH_AMERICA
## 1614             <NA>                            <NA> SOUTH_AMERICA
## 1701             <NA>                            <NA> SOUTH_AMERICA
## 1716             <NA>                            <NA> SOUTH_AMERICA
## 1819 Mammal specimens http://grbio.org/cool/74zs-fm7j NORTH_AMERICA
## 1958      ISM-Mammals                            <NA> NORTH_AMERICA
##      coordinatePrecision coordinateUncertaintyInMeters       country
## 1413                  NA                            NA United States
## 1614                  NA                           505      Colombia
## 1701                  NA                           505      Colombia
## 1716                  NA                           505      Colombia
## 1819                  NA                            20 United States
## 1958                  NA                            NA United States
##      crawlId datasetID                           datasetKey datasetName
## 1413    1454      <NA> 854f602e-f762-11e1-a439-00145eb45e9a        <NA>
## 1614      15      <NA> 87c22676-55df-4e72-b2fa-6488700576d1        <NA>
## 1701      15      <NA> 87c22676-55df-4e72-b2fa-6488700576d1        <NA>
## 1716      15      <NA> 87c22676-55df-4e72-b2fa-6488700576d1        <NA>
## 1819      28      <NA> 06a00852-f764-4fb8-80d4-ca51f0918459        <NA>
## 1958      38      <NA> 07e9980a-a9f6-4a35-9df6-26371459e570        <NA>
##                    dateIdentified day depth depthAccuracy disposition
## 1413                         <NA>   7    NA            NA        <NA>
## 1614                         <NA>  13    NA            NA        <NA>
## 1701                         <NA>  13    NA            NA        <NA>
## 1716                         <NA>  19    NA            NA        <NA>
## 1819 2016-10-24T00:00:00.000+0000  10    NA            NA        <NA>
## 1958                         <NA>  NA    NA            NA        <NA>
##                                dynamicProperties earliestAgeOrLowestStage
## 1413 { "solr_long_lat": "-95.826584,31.230124" }                     <NA>
## 1614                                        <NA>                     <NA>
## 1701                                        <NA>                     <NA>
## 1716                                        <NA>                     <NA>
## 1819                                        <NA>                     <NA>
## 1958                                        <NA>                     <NA>
##      earliestEonOrLowestEonothem earliestEpochOrLowestSeries
## 1413                        <NA>                        <NA>
## 1614                        <NA>                        <NA>
## 1701                        <NA>                        <NA>
## 1716                        <NA>                        <NA>
## 1819                        <NA>                        <NA>
## 1958                        <NA>                        <NA>
##      earliestEraOrLowestErathem earliestPeriodOrLowestSystem elevation
## 1413                       <NA>                         <NA>        NA
## 1614                       <NA>                         <NA>   1827.00
## 1701                       <NA>                         <NA>   1189.00
## 1716                       <NA>                         <NA>   1727.00
## 1819                       <NA>                         <NA>    325.83
## 1958                       <NA>                         <NA>        NA
##      elevationAccuracy endDayOfYear establishmentMeans
## 1413                NA         <NA>               <NA>
## 1614                 0         <NA>               <NA>
## 1701                 0         <NA>               <NA>
## 1716                 0         <NA>               <NA>
## 1819                 0          161               <NA>
## 1958                NA         <NA>             NATIVE
##                         eventDate eventID
## 1413 2015-01-07T00:00:00.000+0000    <NA>
## 1614 2015-09-13T00:00:00.000+0000    <NA>
## 1701 2015-11-13T00:00:00.000+0000    <NA>
## 1716 2015-11-19T00:00:00.000+0000    <NA>
## 1819 2014-06-10T00:00:00.000+0000    <NA>
## 1958 2014-01-01T00:00:00.000+0000    <NA>
##                                                    eventRemarks eventTime
## 1413                                                       <NA>      <NA>
## 1614                                                       <NA>      <NA>
## 1701                                                       <NA>      <NA>
## 1716                                                       <NA>      <NA>
## 1819 OCCURRENCEREMARKS: wild caught; ESTABLISHMENTMEANS: native      <NA>
## 1958                                                       <NA>      <NA>
##           family familyKey fieldNotes fieldNumber formation   fullCountry
## 1413 Dasypodidae      9369       <NA>        0585      <NA> United States
## 1614 Dasypodidae      9369       <NA>        <NA>      <NA>      Colombia
## 1701 Dasypodidae      9369       <NA>        <NA>      <NA>      Colombia
## 1716 Dasypodidae      9369       <NA>        <NA>      <NA>      Colombia
## 1819 Dasypodidae      9369       <NA>        <NA>      <NA> United States
## 1958 Dasypodidae      9369       <NA>        <NA>      <NA> United States
##          gbifID genericName   genus genusKey geodeticDatum
## 1413 1099834716     Dasypus Dasypus  2440775         WGS84
## 1614 1571717883     Dasypus Dasypus  2440775         WGS84
## 1701 1571717875     Dasypus Dasypus  2440775         WGS84
## 1716 1571717891     Dasypus Dasypus  2440775         WGS84
## 1819 1656058026     Dasypus Dasypus  2440775         WGS84
## 1958 1230275499     Dasypus Dasypus  2440775         WGS84
##      geologicalContextID  georeferencedBy     georeferencedDate
## 1413                <NA>             <NA>                  <NA>
## 1614                <NA> Ricardo Ortíz G.            2017-04-26
## 1701                <NA> Ricardo Ortíz G.            2017-04-26
## 1716                <NA> Ricardo Ortíz G.            2017-04-26
## 1819                <NA>          unknown 2016-10-24 00:00:00.0
## 1958                <NA>             <NA>                  <NA>
##                                                                                                                                                                                                                                                                                                                                                                     georeferenceProtocol
## 1413                                                                                                                                                                                                                                                                                                                                                                    digital resource
## 1614 Escobar D, Jojoa LM, Díaz SR, Rudas E, Albarracín RD, Ramírez C, Gómez JY, López CR, Saavedra J, Ortiz R, (2016). Georreferenciación de localidades: Una guía de referencia para colecciones biológicas. Instituto de Investigación de Recursos Biológicos Alexander von Humboldt – Instituto de Ciencias Naturales, Universidad Nacional de Colombia. Bogotá D.C., Colombia. 144 p
## 1701 Escobar D, Jojoa LM, Díaz SR, Rudas E, Albarracín RD, Ramírez C, Gómez JY, López CR, Saavedra J, Ortiz R, (2016). Georreferenciación de localidades: Una guía de referencia para colecciones biológicas. Instituto de Investigación de Recursos Biológicos Alexander von Humboldt – Instituto de Ciencias Naturales, Universidad Nacional de Colombia. Bogotá D.C., Colombia. 144 p
## 1716 Escobar D, Jojoa LM, Díaz SR, Rudas E, Albarracín RD, Ramírez C, Gómez JY, López CR, Saavedra J, Ortiz R, (2016). Georreferenciación de localidades: Una guía de referencia para colecciones biológicas. Instituto de Investigación de Recursos Biológicos Alexander von Humboldt – Instituto de Ciencias Naturales, Universidad Nacional de Colombia. Bogotá D.C., Colombia. 144 p
## 1819                                                                                                                                                                                                                                                                                                                                                                        not recorded
## 1958                                                                                                                                                                                                                                                                                                                                                                                <NA>
##                                                                                                                       georeferenceRemarks
## 1413                                                                                                                                 <NA>
## 1614 Nivel 1. Se valida la coordenada original, coincide con la vereda San Miguel. Incertidumbre por datum y precisión de las coordenadas
## 1701                                    Nivel 1. Se valida la coordenada original. Incertidumbre por datum y precisión de las coordenadas
## 1716                                    Nivel 1. Se valida la coordenada original. Incertidumbre por datum y precisión de las coordenadas
## 1819                                                                                                                                 <NA>
## 1958                                                                                                                                 <NA>
##                                                                                                                                                                                                                                                                                                                                                                                                                                                           georeferenceSources
## 1413                                                                                                                                                                                                                                                                                                                                                                                                                                                             Google Earth
## 1614 Instituto Geográfico Agustín Codazzi - IGAC. Base cartográfica oficial integrada. Plancha 322. Escala 1:100.000. Codificación de la división político administrativa de Colombia(Divipola) - DANE. NASA Land Processes Distributed Active Archive Center (LP DAAC). ASTER Global DEM. EOSDIS/Reverb ECHO. https://reverb.echo.nasa.gov. 2011. Red Nacional de Información RNI, veredas de Colombia, 2016. GeoNames geographical database, http://www.geonames.org/, 2016
## 1701 Instituto Geográfico Agustín Codazzi - IGAC. Base cartográfica oficial integrada. Plancha 281. Escala 1:100.000. Codificación de la división político administrativa de Colombia(Divipola) - DANE. NASA Land Processes Distributed Active Archive Center (LP DAAC). ASTER Global DEM. EOSDIS/Reverb ECHO. https://reverb.echo.nasa.gov. 2011. Red Nacional de Información RNI, veredas de Colombia, 2016. GeoNames geographical database, http://www.geonames.org/, 2016
## 1716 Instituto Geográfico Agustín Codazzi - IGAC. Base cartográfica oficial integrada. Plancha 281. Escala 1:100.000. Codificación de la división político administrativa de Colombia(Divipola) - DANE. NASA Land Processes Distributed Active Archive Center (LP DAAC). ASTER Global DEM. EOSDIS/Reverb ECHO. https://reverb.echo.nasa.gov. 2011. Red Nacional de Información RNI, veredas de Colombia, 2016. GeoNames geographical database, http://www.geonames.org/, 2016
## 1819                                                                                                                                                                                                                                                                                                                                                                                                                                                             not recorded
## 1958                                                                                                                                                                                                                                                                                                                                                                                                                                                                     <NA>
##      georeferenceVerificationStatus habitat
## 1413                           <NA>    <NA>
## 1614                           <NA>    <NA>
## 1701                           <NA>    <NA>
## 1716                           <NA>    <NA>
## 1819                     unverified    <NA>
## 1958          requires verification    <NA>
##                                                                                                                                     higherClassification
## 1413 Animalia; Chordata; Vertebrata; Amniota; Mammalia; Theriiformes-----Theria-Placentalia-Xenarthra; Cingulata; Dasypodoidea; Dasypodidae; Dasypodinae
## 1614                                                                                                                                                <NA>
## 1701                                                                                                                                                <NA>
## 1716                                                                                                                                                <NA>
## 1819                                                                                               Animalia; Chordata; Mammalia; Cingulata; Dasypodidae;
## 1958                                                                                  Animalia | Chordata | Mammalia | Cingulata | Dasypodidae | Dasypus
##                                                           higherGeography
## 1413                  North America; USA; Texas; Leon County; Centerville
## 1614                                                                 <NA>
## 1701                                                                 <NA>
## 1716                                                                 <NA>
## 1819                 North America, United States, Oklahoma, Noble County
## 1958 North America | United States | Missouri | Hickory County |  |  |  |
##      higherGeographyID highestBiostratigraphicZone
## 1413              <NA>                        <NA>
## 1614              <NA>                        <NA>
## 1701              <NA>                        <NA>
## 1716              <NA>                        <NA>
## 1819              <NA>                        <NA>
## 1958              <NA>                        <NA>
##      http://unknown.org/classs http://unknown.org/occurrenceDetails
## 1413                      <NA>                                 <NA>
## 1614                      <NA>                                 <NA>
## 1701                      <NA>                                 <NA>
## 1716                      <NA>                                 <NA>
## 1819                      <NA>                                 <NA>
## 1958                      <NA>                                 <NA>
##      identificationID identificationQualifier identificationReferences
## 1413             <NA>                    <NA>                     <NA>
## 1614             <NA>                    <NA>                     <NA>
## 1701             <NA>                    <NA>                     <NA>
## 1716             <NA>                    <NA>                     <NA>
## 1819             <NA>                       A                     <NA>
## 1958             <NA>                    <NA>                     <NA>
##      identificationRemarks identificationVerificationStatus identifiedBy
## 1413                  <NA>                             <NA>         <NA>
## 1614                  <NA>                             <NA>  L V. García
## 1701                  <NA>                             <NA>     C Guzman
## 1716                  <NA>                             <NA>     C Guzman
## 1819                  <NA>                           legacy      unknown
## 1958                  <NA>                             <NA>         <NA>
##                                                           identifier
## 1413                   urn:uuid:875b3731-d110-4008-bbc3-a21e237e14ce
## 1614                                              UT:CZUT:CZUT-M1764
## 1701                                              UT:CZUT:CZUT-M1766
## 1716                                              UT:CZUT:CZUT-M1770
## 1819 http://arctos.database.museum/guid/UMNH:Mamm:37172?seid=3325000
## 1958                            8b22d40e-7ba4-4d9a-824c-05da299b208e
##      individualCount informationWithheld infraspecificEpithet
## 1413               1                <NA>            mexicanus
## 1614              NA                <NA>                 <NA>
## 1701              NA                <NA>                 <NA>
## 1716              NA                <NA>                 <NA>
## 1819               1                <NA>            mexicanus
## 1958              NA                <NA>                 <NA>
##      institutionCode                                   institutionID
## 1413             YPM                                            <NA>
## 1614              UT                                            <NA>
## 1701              UT                                            <NA>
## 1716              UT                                            <NA>
## 1819            UMNH http://biocol.org/urn:lsid:biocol.org:col:34862
## 1958             ISM http://biocol.org/urn:lsid:biocol.org:col:35004
##      island islandGroup ISO2        key  kingdom kingdomKey language
## 1413   <NA>        <NA>   US 1099834716 Animalia          1       en
## 1614   <NA>        <NA>   CO 1571717883 Animalia          1     <NA>
## 1701   <NA>        <NA>   CO 1571717875 Animalia          1     <NA>
## 1716   <NA>        <NA>   CO 1571717891 Animalia          1     <NA>
## 1819   <NA>        <NA>   US 1656058026 Animalia          1       en
## 1958   <NA>        <NA>   US 1230275499 Animalia          1       en
##                       lastCrawled              lastInterpreted
## 1413 2018-04-23T00:12:16.264+0000 2018-02-03T17:00:27.697+0000
## 1614 2017-11-28T21:02:14.179+0000 2018-02-04T00:48:03.163+0000
## 1701 2017-11-28T21:02:14.179+0000 2018-02-04T00:47:59.235+0000
## 1716 2017-11-28T21:02:14.180+0000 2018-02-04T00:47:58.556+0000
## 1819 2018-05-08T16:58:18.293+0000 2018-02-02T21:18:55.693+0000
## 1958 2017-09-30T14:21:10.054+0000 2018-02-03T21:43:14.977+0000
##                        lastParsed latestEonOrHighestEonothem
## 1413 2017-11-09T20:27:54.444+0000                       <NA>
## 1614 2017-11-28T21:02:14.217+0000                       <NA>
## 1701 2017-11-28T21:02:14.209+0000                       <NA>
## 1716 2017-11-28T21:02:14.199+0000                       <NA>
## 1819 2018-01-09T07:00:05.661+0000                       <NA>
## 1958 2016-12-23T03:21:14.148+0000                       <NA>
##      latestEpochOrHighestSeries latestEraOrHighestErathem
## 1413                       <NA>                      <NA>
## 1614                       <NA>                      <NA>
## 1701                       <NA>                      <NA>
## 1716                       <NA>                      <NA>
## 1819                       <NA>                      <NA>
## 1958                       <NA>                      <NA>
##      latestPeriodOrHighestSystem
## 1413                        <NA>
## 1614                        <NA>
## 1701                        <NA>
## 1716                        <NA>
## 1819                        <NA>
## 1958                        <NA>
##                                                         license lifeStage
## 1413 http://creativecommons.org/publicdomain/zero/1.0/legalcode      <NA>
## 1614       http://creativecommons.org/licenses/by/4.0/legalcode      <NA>
## 1701       http://creativecommons.org/licenses/by/4.0/legalcode      <NA>
## 1716       http://creativecommons.org/licenses/by/4.0/legalcode      <NA>
## 1819 http://creativecommons.org/publicdomain/zero/1.0/legalcode      <NA>
## 1958 http://creativecommons.org/publicdomain/zero/1.0/legalcode      <NA>
##      lithostratigraphicTerms                           locality
## 1413                    <NA>                         Dunn Ranch
## 1614                    <NA>                  Vereda San Miguel
## 1701                    <NA>               Vereda Vega Chiquita
## 1716                    <NA>                 Vereda La Virginia
## 1819                    <NA> 0.6 mi S Perry exit, Interstate 35
## 1958                    <NA>  Route 6, near Jones Spring, Avery
##      locationAccordingTo locationID         locationRemarks
## 1413                <NA>       <NA>                    <NA>
## 1614                <NA>       <NA>                    <NA>
## 1701                <NA>       <NA>                    <NA>
## 1716                <NA>       <NA>                    <NA>
## 1819             unknown       <NA> VERBATIMELEVATION: 1069
## 1958                <NA>       <NA>                    <NA>
##      lowestBiostratigraphicZone member                     modified month
## 1413                       <NA>   <NA>                         <NA>     1
## 1614                       <NA>   <NA>                         <NA>     9
## 1701                       <NA>   <NA>                         <NA>    11
## 1716                       <NA>   <NA>                         <NA>    11
## 1819                       <NA>   <NA> 2017-11-01T14:44:32.000+0000     6
## 1958                       <NA>   <NA> 2014-08-25T11:37:37.000+0000    NA
##      municipality nameAccordingTo namePublishedIn namePublishedInYear
## 1413  Centerville            <NA>            <NA>                <NA>
## 1614         <NA>            <NA>            <NA>                <NA>
## 1701         <NA>            <NA>            <NA>                <NA>
## 1716         <NA>            <NA>            <NA>                <NA>
## 1819         <NA>            <NA>            <NA>                <NA>
## 1958         <NA>            <NA>            <NA>                <NA>
##      nomenclaturalCode
## 1413              ICZN
## 1614              <NA>
## 1701              <NA>
## 1716              <NA>
## 1819              ICZN
## 1958              ICZN
##                                                         occurrenceID
## 1413                   urn:uuid:875b3731-d110-4008-bbc3-a21e237e14ce
## 1614                                              UT:CZUT:CZUT-M1764
## 1701                                              UT:CZUT:CZUT-M1766
## 1716                                              UT:CZUT:CZUT-M1770
## 1819 http://arctos.database.museum/guid/UMNH:Mamm:37172?seid=3325000
## 1958                            8b22d40e-7ba4-4d9a-824c-05da299b208e
##                                                                         occurrenceRemarks
## 1413 MAM number 15982; female female; presonal specimen number DC 20; pregnant, 5 fetuses
## 1614                                                                                 <NA>
## 1701                                                                                 <NA>
## 1716                                                                                 <NA>
## 1819                                                                                 <NA>
## 1958                                                                                 <NA>
##      occurrenceStatus     order orderKey
## 1413             <NA> Cingulata      735
## 1614             <NA> Cingulata      735
## 1701             <NA> Cingulata      735
## 1716             <NA> Cingulata      735
## 1819             <NA> Cingulata      735
## 1958          present Cingulata      735
##                                              organismID organismQuantity
## 1413                                               <NA>             <NA>
## 1614                                               <NA>             <NA>
## 1701                                               <NA>             <NA>
## 1716                                               <NA>             <NA>
## 1819 http://arctos.database.museum/guid/UMNH:Mamm:37172             <NA>
## 1958                                               <NA>             <NA>
##      organismQuantityType organismRemarks originalNameUsage
## 1413                 <NA>            <NA>              <NA>
## 1614                 <NA>            <NA>              <NA>
## 1701                 <NA>            <NA>              <NA>
## 1716                 <NA>            <NA>              <NA>
## 1819                 <NA>            <NA>              <NA>
## 1958                 <NA>            <NA>              <NA>
##            otherCatalogNumbers ownerInstitutionCode parentEventID
## 1413                      <NA>                  YPM          <NA>
## 1614                      <NA>                 <NA>          <NA>
## 1701                      <NA>                 <NA>          <NA>
## 1716                      <NA>                 <NA>          <NA>
## 1819 collector number=EAR 8842                 <NA>          <NA>
## 1958                      <NA>                 <NA>          <NA>
##      parentNameUsage   phylum phylumKey
## 1413            <NA> Chordata        44
## 1614            <NA> Chordata        44
## 1701            <NA> Chordata        44
## 1716            <NA> Chordata        44
## 1819            <NA> Chordata        44
## 1958            <NA> Chordata        44
##                              preparations
## 1413 skeleton (skull only); tissue frozen
## 1614                                 <NA>
## 1701                                 <NA>
## 1716                                 <NA>
## 1819    skin; skull; muscle (95% ethanol)
## 1958                 partial skeleton - 1
##                                                                                      previousIdentifications
## 1413                                                                          Dasypus novemcinctus mexicanus
## 1614                                                                                    Dasypus novemcinctus
## 1701                                                                                    Dasypus novemcinctus
## 1716                                                                                    Dasypus novemcinctus
## 1819 <i>Dasypus novemcinctus mexicanus</i> (accepted ID) identified by unknown on 2016-10-24; method: legacy
## 1958                                                                                                    <NA>
##         protocol publishingCountry                     publishingOrgKey
## 1413 DWC_ARCHIVE                US 2e167bb0-4441-11db-9ba2-b8a03c50a862
## 1614 DWC_ARCHIVE                CO 5a45153b-bdf9-44ae-b7a7-e3261896540b
## 1701 DWC_ARCHIVE                CO 5a45153b-bdf9-44ae-b7a7-e3261896540b
## 1716 DWC_ARCHIVE                CO 5a45153b-bdf9-44ae-b7a7-e3261896540b
## 1819 DWC_ARCHIVE                US 2d6267a0-0561-11d8-b851-b8a03c50a862
## 1958 DWC_ARCHIVE                US b0e368b3-da4b-45ff-8d24-12ffd422bc38
##                                                       recordedBy
## 1413 Gunter P. Wagner, Wess Dunn, Kendell Dunn, Mihaela Pavlicev
## 1614                                                 L V. García
## 1701                                                    C Guzman
## 1716                                                    C Guzman
## 1819                               Collector(s): Eric A. Rickart
## 1958                                                Widga, Henry
##      recordNumber
## 1413         <NA>
## 1614         <NA>
## 1701         <NA>
## 1716         <NA>
## 1819     EAR 8842
## 1958         <NA>
##                                                            references
## 1413 http://collections.peabody.yale.edu/search/Record/YPM-MAM-015982
## 1614                                                             <NA>
## 1701                                                             <NA>
## 1716                                                             <NA>
## 1819               http://arctos.database.museum/guid/UMNH:Mamm:37172
## 1958            http://portal.vertnet.org/o/ism/ism-mammals?id=693990
##      reproductiveCondition rights                           rightsHolder
## 1413   pregnant, 5 fetuses   <NA> Yale Peabody Museum of Natural History
## 1614                  <NA>   <NA>                                   <NA>
## 1701                  <NA>   <NA>                                   <NA>
## 1716                  <NA>   <NA>                                   <NA>
## 1819                  <NA>   <NA>                                   <NA>
## 1958                  <NA>   <NA>                                   <NA>
##      samplingEffort samplingProtocol
## 1413           <NA>             <NA>
## 1614           <NA>     Huellla Yeso
## 1701           <NA>     Huellla Yeso
## 1716           <NA>     Huellla Yeso
## 1819           <NA>             <NA>
## 1958           <NA>             <NA>
##                                   scientificName scientificNameID  sex
## 1413 Dasypus novemcinctus mexicanus Peters, 1864             <NA> <NA>
## 1614         Dasypus novemcinctus Linnaeus, 1758             <NA> <NA>
## 1701         Dasypus novemcinctus Linnaeus, 1758             <NA> <NA>
## 1716         Dasypus novemcinctus Linnaeus, 1758             <NA> <NA>
## 1819 Dasypus novemcinctus mexicanus Peters, 1864             <NA> <NA>
## 1958         Dasypus novemcinctus Linnaeus, 1758             <NA> <NA>
##                   species speciesKey specificEpithet startDayOfYear
## 1413 Dasypus novemcinctus    2440779    novemcinctus           <NA>
## 1614 Dasypus novemcinctus    2440779    novemcinctus           <NA>
## 1701 Dasypus novemcinctus    2440779    novemcinctus           <NA>
## 1716 Dasypus novemcinctus    2440779    novemcinctus           <NA>
## 1819 Dasypus novemcinctus    2440779    novemcinctus           <NA>
## 1958 Dasypus novemcinctus    2440779    novemcinctus           <NA>
##      taxonConceptID taxonID taxonKey taxonomicStatus  taxonRank
## 1413           <NA>    <NA>  6163121            <NA> SUBSPECIES
## 1614           <NA>    <NA>  2440779        accepted    SPECIES
## 1701           <NA>    <NA>  2440779        accepted    SPECIES
## 1716           <NA>    <NA>  2440779        accepted    SPECIES
## 1819           <NA>    <NA>  6163121            <NA> SUBSPECIES
## 1958           <NA>    <NA>  2440779            <NA>    SPECIES
##                                   taxonRemarks           type typeStatus
## 1413 Animals and Plants: Vertebrates - Mammals PhysicalObject       <NA>
## 1614                                      <NA>  Objeto físico       <NA>
## 1701                                      <NA>  Objeto físico       <NA>
## 1716                                      <NA>  Objeto físico       <NA>
## 1819                                      <NA> PhysicalObject       <NA>
## 1958                                      <NA> PhysicalObject       <NA>
##      typifiedName verbatimCoordinateSystem verbatimElevation
## 1413         <NA>                     <NA>              <NA>
## 1614         <NA>                     <NA>              1871
## 1701         <NA>                     <NA>               950
## 1716         <NA>                     <NA>              1473
## 1819         <NA>          decimal degrees              <NA>
## 1958         <NA>                     <NA>              <NA>
##      verbatimEventDate
## 1413              <NA>
## 1614              <NA>
## 1701              <NA>
## 1716              <NA>
## 1819         6/10/2014
## 1958     July 22, 2014
##                                               verbatimLocality verbatimSRS
## 1413                                                      <NA>        <NA>
## 1614                                                San Miguel        <NA>
## 1701                                             Vega chiquita        <NA>
## 1716                                               La Virginia        <NA>
## 1819                          United States | Oklahoma | Noble        <NA>
## 1958 North America | United States | Missouri | Hickory County        <NA>
##      verbatimTaxonRank
## 1413              <NA>
## 1614              <NA>
## 1701              <NA>
## 1716              <NA>
## 1819              <NA>
## 1958              <NA>
##                                                                                          vernacularName
## 1413 Nine-banded Armadillo; long-nosed armadillos; Armadillos; mammals; vertebrates; chordates; animals
## 1614                                                                                               <NA>
## 1701                                                                                               <NA>
## 1716                                                                                               <NA>
## 1819                                                                                               <NA>
## 1958                                                                                               <NA>
##      waterBody year downloadDate cut_block
## 1413      <NA> 2015   2018-05-15         4
## 1614      <NA> 2015   2018-05-15         2
## 1701      <NA> 2015   2018-05-15         2
## 1716      <NA> 2015   2018-05-15         2
## 1819      <NA> 2014   2018-05-15         3
## 1958      <NA> 2014   2018-05-15         4
```

```r
#occ_train <- subset(occ_final,cut_block!=4) # this is the selection to be used for model training
#occ_test <- subset(occ_final,cut_block==4) # this is the opposite of the selection which will be used for model testing
plot(occ_final)
plot(subset(occ_final,cut_block==1),col=1,add=T)
plot(subset(occ_final,cut_block==2),col=2,add=T)
plot(subset(occ_final,cut_block==3),col=3,add=T)
plot(subset(occ_final,cut_block==4),col=4,add=T)
```

![](Appendix1_case_study_files/figure-html/cut_occ_into_training_testing-2.png)<!-- -->

```r
#cut_random <- get.checkerboard2(occ=as.data.frame(occ@coords), bg.coords=as.data.frame(bg@coords), aggregation.factor = c(2, 2), env=env)
#occ_update@data$cut_random <- cut_random$occ.grp
#bg@data$cut_random <- cut_random$bg.grp
```

###2.5 Format data for Maxent
The data input can either be spatial or tabular. In our example, we use the tabular format, which can be potentially more flexible. We extract environmental conditions for background, training, and testing points in a dataframe format. 

#####Thread 15

```r
# extracting env conditions for training occ from the raster stack;
# a data frame is returned (i.e multiple columns)
p <- extract(clim,occ_train) 
# env conditions for testing occ
p_test <- extract(clim,occ_test) 
# extracting env conditions for background
a <- extract(clim,bg)  
```

Maxent reads a "1" as presence and "0" as pseudo-absence. Thus, we need to assign a "1" to the training environmental conditions and a "0" for the background. We create a set of rows with the same number as the training and testing data, and put the value of "1" for each cell and a "0" for background. We combine the "1"s and "0"s into a vector that was added to the dataframe containing the environmental conditions associated with the testing and background conditions.

#####Thread 16

```r
# repeat the number 1 as many numbers as the number of rows in p, 
# and repeat 0 as the rows of background points
pa <- c(rep(1,nrow(p)), rep(0,nrow(a))) 

# (rep(1,nrow(p)) creating the number of rows as the p data set to 
# have the number '1' as the indicator for presence; rep(0,nrow(a)) 
# creating the number of rows as the a data set to have the number
# '0' as the indicator for absence; the c combines these ones and 
# zeros into a new vector that can be added to the Maxent table data
# frame with the environmental attributes of the presence and absence locations
pder <- as.data.frame(rbind(p,a)) 
```

##3 Maxent models
###3.1 Simple implementation

#####Thread 17

```r
# train Maxent with spatial data
# mod <- maxent(x=clim,p=occ_train)

# train Maxent with tabular data
mod <- dismo::maxent(x=pder, ## env conditions
              p=pa,   ## 1:presence or 0:absence
              path=paste0(getwd(),"/output/maxent_outputs"), ## folder for maxent output; 
              # if we do not specify a folder R will put the results in a temp file, 
              # and it gets messy to read those. . .
              args=c("responsecurves") ## parameter specification
              )
# the maxent functions runs a model in the default settings. To change these parameters,
# you have to tell it what you want...i.e. response curves or the type of features

# view the maxent model in a html brower
mod
```

```
## class    : MaxEnt 
## variables: bio1 bio10 bio11 bio12 bio13 bio14 bio15 bio16 bio17 bio18 bio19 bio2 bio3 bio4 bio5 bio6 bio7 bio8 bio9
```

```r
# view detailed results
# mod@results
```

###3.2 Predict function
Running Maxent in R will not automatically make projection to layers, unless you specify this using the parameter *projectionlayers*. However, we could make projections (to dataframes or raster layers) post hoc using the *predict* function.

#####Thread 18

```r
# example 1, project to study area [raster]
ped1 <- predict(mod,studyArea) # studyArea is the clipped rasters 
plot(ped1) # plot the continuous prediction
```

![](Appendix1_case_study_files/figure-html/predict-1.png)<!-- -->

```r
# example 2, project to the world
#ped2 <- predict(mod,clim)
#plot(ped2)

# example 3, project with training occurrences [dataframes]
ped3 <- predict(mod,p)
head(ped3)
```

```
## [1] 0.8311057 0.3677974 0.8376744 0.8308824 0.3057640 0.8241310
```

```r
hist(ped3)# creates a histogram of the prediction
```

![](Appendix1_case_study_files/figure-html/predict-2.png)<!-- -->

###3.3 Model evaluation
To evaluate models, we use the *evaluate* function from the "dismo" package. Evaluation indices include AUC, TSS, Sensitivity, Specificity, etc.

#####Thread 19

```r
# using "training data" to evaluate 
#p & a are dataframe/s (the p and a are the training presence and background points)
mod_eval_train <- dismo::evaluate(p=p,a=a,model=mod) 
print(mod_eval_train)
```

```
## class          : ModelEvaluation 
## n presences    : 267 
## n absences     : 10000 
## AUC            : 0.8601406 
## cor            : 0.2621316 
## max TPR+TNR at : 0.4345195
```

```r
mod_eval_test <- dismo::evaluate(p=p_test,a=a,model=mod)  
print(mod_eval_test) # training AUC may be higher than testing AUC
```

```
## class          : ModelEvaluation 
## n presences    : 267 
## n absences     : 10000 
## AUC            : 0.7974086 
## cor            : 0.2135319 
## max TPR+TNR at : 0.4066744
```

To threshold our continuous predictions of suitability into binary predictions we use the threshold function of the "dismo" package. To plot the binary prediction, we plot the predictions that are larger than the threshold.  

#####Thread 20

```r
# calculate thresholds of models
thd1 <- threshold(mod_eval_train,"no_omission") # 0% omission rate 
thd2 <- threshold(mod_eval_train,"spec_sens") # highest TSS

# plotting points that are above the previously calculated thresholded value
plot(ped1>=thd1) 
```

![](Appendix1_case_study_files/figure-html/model_evaluation2-1.png)<!-- -->

##4 Maxent parameters
###4.1 Select features

#####Thread 21

```r
# load the function that prepares parameters for maxent
source("Appendix2_prepPara.R")

mod1_autofeature <- maxent(x=pder[c("bio1","bio4","bio11")], 
                           ## env conditions, here we selected only 3 predictors
                           p=pa,
                           path=paste0(getwd(),"/output/maxent_outputs"),
                           ## 1:presence or 0:absence
                           #path=, # use preppara to prepare the path
                           ## this is the folder you will find manxent output
                           args=prepPara(userfeatures=NULL) 
                           ) 
                           ## default is autofeature
              
# or select Linear& Quadratic features
mod1_lq <- maxent(x=pder[c("bio1","bio4","bio11")],
                  p=pa,
                  path=paste0(getwd(),"/output/maxent_outputs1_lq"),
                  args=prepPara(userfeatures="LQ") ) 
                  ## default is autofeature, here LQ represents Linear& Quadratic
                  ## (L-linear, Q-Quadratic, H-Hinge, P-Product, T-Threshold)
```

###4.2 Change beta-multiplier

#####Thread 22

```r
#change betamultiplier for all features
mod2 <- maxent(x=pder[c("bio1","bio4","bio11")], 
               p=pa, 
              path=paste0(getwd(),"/output/maxent_outputs2_0.5"), 
              args=prepPara(userfeatures="LQ",
                            betamultiplier=0.5) ) 

mod2 <- maxent(x=pder[c("bio1","bio4","bio11")], 
               p=pa, 
              path=paste0(getwd(),"/output/maxent_outputs2_complex"), 
              args=prepPara(userfeatures="LQH",
                            ## include L, Q, H features
                            beta_lqp=1.5, 
                            ## use different betamultiplier for different features
                            beta_hinge=0.5 ) ) 
```

###4.3 Specify projection layers

#####Thread 23

```r
# note: (1) the projection layers must exist in the hard disk (as relative to computer RAM); 
# (2) the names of the layers (excluding the name extension) must match the names 
# of the predictor variables; 
mod3 <- maxent(x=pder[c("bio1","bio11")], 
               p=pa, 
              path=paste0(getwd(),"/output/maxent_outputs3_prj1"), 
              args=prepPara(userfeatures="LQ",
                            betamultiplier=1,
                            projectionlayers="/data/studyarea") ) 

# load the projected map
ped <- raster(paste0("output/maxent_outputs3_prj1/species_studyarea.asc"))
plot(ped)

# we can also project on a broader map, but please 
# caustion about the inaccuracy associated with model extrapolation.
mod3 <- maxent(x=pder[c("bio1","bio11")], 
               p=pa, 
              path=paste0(getwd(),"/output/maxent_outputs3_prj2"), 
              args=prepPara(userfeatures="LQ",
                            betamultiplier=1,
                            projectionlayers="/data/studyarea") ) 
# plot the map
ped <- raster(paste0("output/maxent_outputs3_prj2/species_studyarea.asc"))
plot(ped)
```

![](Appendix1_case_study_files/figure-html/specify_projection_layers_data_table-1.png)<!-- -->

```r
# simply check the difference if we used a different betamultiplier
mod3_beta1 <- maxent(x=pder[c("bio1","bio11")], 
               p=pa, 
              path=paste0(getwd(),"/output/maxent_outputs3_prj3"), 
              args=prepPara(userfeatures="LQ",
                            betamultiplier=100, 
                            ## for an extreme example, set beta as 100
                            projectionlayers="/data/bioclim") ) 
```

```
## Warning in matrix(as.numeric(d)): NAs introduced by coercion
```

```r
ped3 <- raster(paste0("output/maxent_outputs3_prj3/species_bioclim.asc"))

plot(ped-crop(ped3,extent(ped))) ## quickly check the difference between the two predictions
```

![](Appendix1_case_study_files/figure-html/specify_projection_layers_data_table-2.png)<!-- -->

###4.4 Clamping function

#####Thread 24

```r
# enable or disable clamping function; note that clamping function is involved when projecting
mod4_clamp <- maxent(x=pder[c("bio1","bio11")],
                     p=pa,
                     path=paste0(getwd(),"/output/maxent_outputs4_clamp"), 
                     args=prepPara(userfeatures="LQ",
                                   betamultiplier=1,
                                   doclamp = TRUE,
                                   projectionlayers="/data/bioclim")) 

mod4_noclamp <- maxent(x=pder[c("bio1","bio11")], 
                       p=pa, 
                       path=paste0(getwd(),"/output/maxent_outputs4_noclamp"),
                       args=prepPara(userfeatures="LQ",
                                      betamultiplier=1,
                                      doclamp = FALSE,
                                      projectionlayers="/data/bioclim") ) 

ped_clamp <- raster(paste0("output/maxent_outputs4_clamp/species_bioclim.asc") )
ped_noclamp <- raster(paste0("output/maxent_outputs4_noclamp/species_bioclim.asc") )
plot(stack(ped_clamp,ped_noclamp))
```

![](Appendix1_case_study_files/figure-html/clamping_function-1.png)<!-- -->

```r
plot(ped_clamp - ped_noclamp) 
```

![](Appendix1_case_study_files/figure-html/clamping_function-2.png)<!-- -->

```r
## we may notice small differences, especially clamp shows higher predictions in most areas.
```

###4.5 Cross validation

#####Thread 25

```r
mod4_cross <- maxent(x=pder[c("bio1","bio11")], p=pa, 
                            path=paste0(getwd(),"/output/maxent_outputs4_cross"), 
                            args=prepPara(userfeatures="LQ",
                                          betamultiplier=1,
                                          doclamp = TRUE,
                                          projectionlayers="/data/bioclim",
                                          replicates=5, ## 5 replicates
                                          replicatetype="crossvalidate") )
                                          ##possible values are: crossvalidate,bootstrap,subsample
```

##5 experiment of "high AUC" model

```r
# experiment 1: smaller training area vs. whole continent

clim_list <- list.files("data/bioclim/",pattern=".bil$",full.names = T) # '..' leads to the path above the folder where the .rmd file is located
# stacking the bioclim variables to process them at one go 
clim <- raster::stack(clim_list) 

# model for smaller training area # model 1: smaller training area, trying to simulate accessible area (M)
small_extent <- buffer(occ_final,width=100000) #unit is meter 

#country_shp <- shapefile("data/GIS polygon/country.shp")
#extent(occ_final)
big_extent <- extent(c(-130,-20,-60,60))

small_area <- mask(clim,  small_extent)
big_area   <- crop(clim,  big_extent)
plot(small_area[[1]])
```

![](Appendix1_case_study_files/figure-html/experiment-1.png)<!-- -->

```r
plot(big_area[[1]])
```

![](Appendix1_case_study_files/figure-html/experiment-2.png)<!-- -->

```r
dir.create("data/big_extent")
```

```
## Warning in dir.create("data/big_extent"): 'data/big_extent' already exists
```

```r
dir.create("data/small_extent")
```

```
## Warning in dir.create("data/small_extent"): 'data/small_extent' already
## exists
```

```r
writeRaster(small_area,
            # a series of names for output files
            filename=paste0("data/small_extent/",names(small_area),".asc"), 
            format="ascii", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)
writeRaster(big_area,
            # a series of names for output files
            filename=paste0("data/big_extent/",names(big_area),".asc"), 
            format="ascii", ## the output format
            bylayer=TRUE, ## this will save a series of layers
            overwrite=T)

set.seed(1) 
small_bg <- sampleRandom(x=small_area,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 
set.seed(1) 
big_bg <- sampleRandom(x=big_area,
                   size=10000,
                   na.rm=T, #removes the 'Not Applicable' points  
                   sp=T) # return spatial points 

p <- extract(clim,occ_train) 
a_small <- extract(clim,small_bg)  
a_big <- extract(clim,big_bg)  

pa_small <- c(rep(1,nrow(p)), rep(0,nrow(a_small))) 
pa_big   <- c(rep(1,nrow(p)), rep(0,nrow(a_big))) 

pder_small <- as.data.frame(rbind(p,a_small)) 
pder_big   <- as.data.frame(rbind(p,a_big)) 


small_mod <- maxent(x=pder_small,#[c("bio1","bio11")], 
                    p=pa_small, 
              path=paste0(getwd(),"/output/experiment_small"), 
              args=prepPara(projectionlayers="/data/big_extent") 
              )
big_mod <- maxent(x=pder_big,#[c("bio1","bio11")], 
                  p=pa_big, 
                  path=paste0(getwd(),"/output/experiment_big"), 
                  args=prepPara(projectionlayers="/data/big_extent") )

# show the diff of small/big model
ped_small <- raster("output/experiment_small/species_big_extent.asc")
ped_big   <- raster("output/experiment_big/species_big_extent.asc")
ped_combined <- stack(ped_small,ped_big)
names(ped_combined) <- c("small_trainingArea","big_trainArea")
plot( ped_combined )
```

![](Appendix1_case_study_files/figure-html/experiment-3.png)<!-- -->

```r
# compare the AUC , extract test occ data, 
p_test <- extract(clim,occ_test) 

mod_eval_small <- dismo::evaluate(p=p_test,a=a_small,model=small_mod) 
print(mod_eval_small)
```

```
## class          : ModelEvaluation 
## n presences    : 267 
## n absences     : 10000 
## AUC            : 0.6528989 
## cor            : 0.09312594 
## max TPR+TNR at : 0.472774
```

```r
mod_eval_big   <- dismo::evaluate(p=p_test,a=a_big,  model=big_mod) 
print(mod_eval_big)
```

```
## class          : ModelEvaluation 
## n presences    : 267 
## n absences     : 10000 
## AUC            : 0.9054127 
## cor            : 0.3391529 
## max TPR+TNR at : 0.3146109
```








