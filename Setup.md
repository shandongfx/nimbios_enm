# Overview

This tutorial will guide you to process spatial data and implement Maxent modeling in R. Therefore, you need to make sure necessary software/libraries are ready in your computer.

### Step 1: install the list of software and libraries 
(partly adapted from [data carpentry](https://datacarpentry.org/geospatial-workshop/setup.html))

| Name | Type | Need | Install | Description and notes |
| -------- | ------------ | - | ------------- | ----------- |
| R | software | Yes | Go to CRAN’s [cloud download page](https://cloud.r-project.org/) and select the version for your operating system. Download the base subdirectory and install it on your computer. | Software environment for statistical and scientific computing |
| RStudio | software | Yes| Go to RStudio [download page](https://www.rstudio.com/products/rstudio/download/#download) and select RStudio Desktop for your operating system. Download and install it on your computer.| Graphic User Interface (GUI) for R|
|Java JDK| software | Yes | Go to JAVA [download page](https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) and select Java SE Development Kit 8u191 for your operating system.  Download and install version x86 if your operating system is 32-bit, or use version x64 if your operating system is 64-bit. Read more about [32-bit vs. 64-bit](https://en.wikipedia.org/wiki/X86-64). |A programing environmental to run [Maxent](https://github.com/mrmaxent/Maxent/tree/master/ArchivedReleases) modeling algorithm. This is different than Java (which you might already have installed on your computer)|
|raster|R package|Yes|Option 1: install packages from Rstudio interface.<br><br>Option 2: use install.packages() function in R terminal.  <br><br>Option 3: run the following script in R terminal to install:<br> <code>packages_needed <- c("raster","rgdal","sp","dismo","ENMeval", "jsonlite" , "rJava")  <br> pk_to_install <- packages_needed [!( packages_needed %in% rownames(installed.packages())  )]  <br>  if( length(pk_to_install) > 0 ) {  <br>install.packages(pk_to_install,repos="http://cran.r-project.org")  <br> }  </code>  |for raster analysis|
|sp |R package|Yes|see above|for spatial analysis|
|rgdal|R package|Yes|see above|for spatial analysis|
|dismo|R package|Yes|see above|a collection of ENM/SDM tools, including a function to run Maxent.jar in R|
|rJava|R package|Yes|Note: in macOS, if you see error like this when loading rJava library:<br>error: unable to load shared object '/Library/Frameworks/R.framework/Versions/3.5/Resources/library/rJava/libs/rJava.so':  Please run the following code in your terminal:  <code>sudo R CMD javareconf</code>|An interface to Java|
|jsonlite|R package|Yes|see above|necessary for download data from GBIF|
|ENMeval|R package|Yes|see above|a collection of ENM/SDM tools, including a function to separate occurrences|
|Maxent|software|Yes|Install dismo package first, then<br>Option 1: manually download from this [link](https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar) and move Maxent.jar to the path where dismo package is installed, which can be obtained from this function<code>system.file("java", package="dismo")</code>.  <br><br>Option 2: run the following script in R terminal to download.<br><code>if( !file.exists(paste0(system.file("java", package="dismo"),"/maxent.jar"))  )   {<br>utils::download.file(url="https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar",<br>destfile=paste0(system.file("java", package="dismo"),"/maxent.jar"),<br>mode="wb")<br>}</code>|Maxent modeling algorithm|
|GDAL|software|optional|Windows: do the installations through OSGeo4W. <br><code>-Download either the [32 bit or 64 bit OSGeo4W installer](https://trac.osgeo.org/osgeo4w/). -Run the OSGeo4W setup program.<br>-Select “Advanced Install” and press Next.<br>-Select “Install from Internet” and press Next.<br>-Select a installation directory. The default suggestion is fine in most cases. Press Next.<br>-Select “Local packacke directory”. The suggestions is fine in most cases. Press Next.<br>-Select “Direct connection” and press Next.<br>-Choose the download.osgeo.org and press Next.<br>-Find “gdal” or “proj” under “Commandline_Utilities” and click the package in the “New” column until the version you want to install appears.<br>-Press next to install PROJ.</code> <br><br>macOS: run the following code.<br><code>$ brew tap osgeo/osgeo4mac && brew tap --repair<br>$ brew install proj<br>$ brew install geos<br>$ brew install gdal2 --with-armadillo --with-complete --with-libkml --with-unsupported<br>$ brew link --force gdal2</code><br><br>Find more help via this [link](https://proj4.org/install.html).| Geospatial model for reading and writing a variety of formats; this is necessary if you want to install rgdal package from source code|
|PROJ.4|software|optional|see above|Coordinate reference system transformations; this is necessary if you want to install rgdal package from source code|

   
   
### Step 2: test if you have successfully install all necessary software/packages:
#### 2.1 Open RStudio
#### 2.2 Load libraries from R terminal:
```R
library("raster")
library("dismo")
library("rgdal")
library("sp")
library("ENMeval")
library("rJava")
```
#### 2.3 Run a simple Maxent model from R terminal:
```R
# get predictor variables
fnames <- list.files(path=paste(system.file(package="dismo"), '/ex', sep=''), 
              pattern='grd', full.names=TRUE )
predictors <- stack(fnames)

# file with presence points
occurence <- paste(system.file(package="dismo"), '/ex/bradypus.csv', sep='')
occ <- read.table(occurence, header=TRUE, sep=',')[,-1]

# witholding a 20% sample for testing 
fold <- kfold(occ, k=5)
occtest <- occ[fold == 1, ]
occtrain <- occ[fold != 1, ]

# fit model, biome is a categorical variable
me <- maxent(predictors, occtrain, factors='biome')

# see the maxent results in a browser:
me
```

