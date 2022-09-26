# MSc Dissertation
This is a repository which contains some of the main files written during my MSc Dissertation, 'Benchmarking predictions of post-disaster human displacement'.

## ODDRIN R Files
Please note that some of the code presented here depends on the installation of ODDRIN, the downloading and storing of relevant data, and some packages and functions in the ODDRIN files. These functions are mostly not included in the files given here unless absolutely necessary, with the aim of distinguishing between code written during this dissertation and that which was already created. The ODDRIN repository, which contains installation instructions and some files which are referenced below, [can be found here.]

Roughly speaking, the first three files below (ODDpolys, DataViz, and DataPrep) are new, the following two files (DispXNew and BDXNew) are new versions of previously created files which have been developed during this dissertation, and the last two files (ModelChanges and RunAlgo) contain only minor edits to previously created functions and files and are included here for completeness. Some of the new files contain functions which build off previous versions of functions, and we do our best to note whenever this is the case below. 

![plot](./Sanma.jpeg)

##### ODDpolys.R

Key Functions:
- `convIso3Country` - function was already built in ODDRIN which converts a country's ISO3 code to the country name.
- `cleanValData` - filters to the inputted country iso3 code, removes NA columns, orders the rows such that highest admin levels are at the  bottom, formats the date columns, adds a column describing the polygon which will be used in getPolys to extract data from OSM.
- `getPolys` - extracts polygons from OSM given the validation data spreadsheet.
- `extractIndices` - extracts the indices of the rows of the coordinates in the ODDpixels object that have points contained within the polygons that were extracted from getPolys.
- `initialize ODDpolys` - fill the slots of ODDpolys with polyIndices (a list of vectors, one for each polygon, with each vector containing the indices of the ODDpixels object whose points lie within the polygon), sourceInfo (data frame made up of source date, source type, and the source itself), valDF (dataframe with the validation data, i.e. mortality, displacement, buildings damaged, buildings destroyed), polygons (the polygons extracted from OSM), bbox (the bounding box of the polygons), and data (the data slot of the polys Spatial Polygons Data Frame extracted from OSM using getPolys).
- `inPoly`: function which was already built in ODDRIN which allows the extraction of point indices within a grid.
- `ExtractOSMBuild` - function which was already built in ODDRIN which allows extraction of OSM buildings for a given bbox. Warning: this function may take several minutes to run, depending on the size of the bbox.
- `ParAggnBuildings` - function which was already built in ODDRIN which allows the aggregation of points to pixels. Warning: this function will take minutes or hours depending on the size of the bbox.

##### DataViz.R

This file contains much of the code that was used to produce the plots throughout the dissertation. 

##### DataPrep.R

This file contains the inference of unaffected buildings discussed in Section 4.2.1 of the dissertation, along with the interpolation of population and hazard intensity values. Note that downloading old OpenStreetMap from before Cyclone Harold is required; this data can be found above and is entitled 'hotosm_vut_buildings.gpkg'.

##### DispXNew.R

Key Functions:
- `DispX_new` - a new version of DispX - a function which was already built in ODDRIN which predicts the number of people displaced, deceased, and buildings destroyed conditional on the hazard, model, and parameterisation. The main changes for this new version are the aggregation of predictions at a polygonal level (this aggregation depends on the ODDpolys.R file and data) and the extension to predict the quantity of buildings damaged and unaffected.
- `LL_IDP_new` - function which presents the output of the kernel density function, which was already built in ODDRIN, but has been extended slightly to work with polygons and with buildings damaged.

##### BDXNew.R

Key Functions:
- `BDX_new` - a new version of BDX - a function which was already built in ODDRIN which classifies the buildings in the Copernicus data according to different building damage gradings. This version of BDX uses the MLR method described in Section 4.2.3 in the dissertation. 

We also include an example at the end of this file which shows how we conducted the 5-fold cross-validation. This example requires using BDX_new to obtain the UnscaledVals for the given Omega parameterisation first, before being able to run the 5-fold cross-validation.

##### ModelChanges.R

This file is included for completeness and only contains minor changes from functions already present in ODDRIN. It is important to include as it shows the calculation of the objective function (i.e. the distance for the ABC-SMC algorithm).

##### RunAlgo.R

Much the same as above; this file is included for completeness and contains only minor modifications to some functions which were already present in ODDRIN. Running only this file as well as those above allows the running of the ABC-SMC algorithm.

## CLIMADA Python Files

![plot](./IntensityMap.png)

Similar to the above, running this CLIMADA file will require the installation of both Python and CLIMADA's software, instructions for which [are available here.]

##### HaroldCLIMADA.py

This file contains the code which is used for much of the visualisations that are cited as using CLIMADA's functionality throughout the dissertation. The code also steps through exposure and impact calculations using CLIMADA's model for Cyclone Harold. Note that CLIMADA uses an economic cost approach for these calculations.

## Data

The Cyclone Harold data can be found in tc_harold_val_data_v6.csv. Some of the main objects have been provided on this repository (excluding those which were too large in size), as well as the old OpenStreetMap data. Other relevant data can be found on the ODDRIN repository, as well as through citations in Section 2.2 of the dissertation.

## Comments
Some of the code and functions can take several days or even weeks to run, depending on parameter choices. When this is the case, we have tried to explicitly warn the user and comment out very intense computations. The full ODDpolys and BD objects which are used in calculations are provided, to save the user time in initialising these. The ODDpixels object was too large to store here. For data requests or queries, please email james-obrien1@outlook.com or contact [@JamesOBrien.]

[can be found here.]: https://github.com/hamishwp/ODDRIN

[are available here.]: https://climada-python.readthedocs.io/en/v3.2.0/guide/Guide_Installation.html

[@JamesOBrien.]: https://www.linkedin.com/in/james-obrien1/
