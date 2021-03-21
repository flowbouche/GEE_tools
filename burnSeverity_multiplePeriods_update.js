/*

Burn Severity Mapping: multitemporal analysis of burn severity and detection of burned areas

Adam Bouché

This script provides a function allowing for burned severity analysis for multiple time periods, modifying and
expanding on a script originally published by the UN-SPIDER program.

A spectral index called the Normalized Burn Ratio is calculated for the pre-fire and post-fire satellite images
for each analysis period using the Near Infrared (NIR) and Short-wave infrared (SWIR2) bands. Raster algebra is
then performed to calculate the change in the NBR values (dNBR). This script is designed for processing images
from Landsat-8 and Sentinel-2 sensors.

The single-analysis-period script on which this script is based can be found here:
https://un-spider.org/advisory-support/recommended-practices/recommended-practice-burn-severity/burn-severity-earth-engine

The source for the method is:
Petropoulos, G. P., Griffiths, H. M., & Kalivas, D. P. (2014). Quantifying spatial and temporal vegetation
recovery dynamics following a wildfire event in a Mediterranean landscape using EO data and GIS.
Applied Geography, 50, 120–131.
*/


//==========================================================================================
//   ANALYSIS SETTINGS - ATTN USER!!! input analysis settings here
//==========================================================================================

// INITIALIZE VARIABLES USED IN THE FUNCTIONS

// USER: CHOOSE SATELLITE IMAGERY SOURCE
// var platform = 'L8' // For Landsat-8
var platform = 'S2' // For Sentinel-2

// USER: DEFINE AN AREA OF INTEREST (AOI)
var aoi = ee.Geometry.Polygon([
  [
    [ -72.435883, -35.540058 ], [ -72.431499, -35.544763 ], [ -72.426375, -35.544077 ], [ -72.424671, -35.539961 ], [ -72.423897, -35.533738 ], [ -72.427985, -35.534038 ], [ -72.431896, -35.530179 ], [ -72.437491, -35.53002  ], [ -72.440437, -35.527437 ], [ -72.444927, -35.525227 ], [ -72.444293, -35.522331 ], [ -72.433256, -35.514318 ], [ -72.424906, -35.509559 ], [ -72.414699, -35.509016 ], [ -72.407876, -35.504213 ], [ -72.402089, -35.499796 ], [ -72.394356, -35.497516 ], [ -72.389394, -35.488496 ], [ -72.382101, -35.484537 ], [ -72.376473, -35.483862 ], [ -72.370793, -35.481939 ], [ -72.365708, -35.482081 ], [ -72.353893, -35.479495 ], [ -72.344741, -35.47975  ], [ -72.341638, -35.478587 ], [ -72.332675, -35.483415 ], [ -72.320215, -35.477513 ], [ -72.316775, -35.480522 ], [ -72.306029, -35.479152 ], [ -72.301453, -35.479277 ], [ -72.293792, -35.478654 ], [ -72.251728, -35.508528 ], [ -72.248369, -35.513615 ], [ -72.247737, -35.523209 ], [ -72.248498, -35.529435 ], [ -72.246188, -35.535326 ], [ -72.245506, -35.543673 ], [ -72.243163, -35.548733 ], [ -72.231825, -35.558198 ], [ -72.23094,  -35.561553 ], [ -72.226918, -35.56291  ], [ -72.227561, -35.566224 ], [ -72.223113, -35.569674 ], [ -72.22323, -35.572586  ], [ -72.22343,  -35.577577 ], [ -72.226209, -35.583332 ], [ -72.221184, -35.585133 ], [ -72.214532, -35.584478 ],
    [ -72.209506, -35.586277 ], [ -72.209606, -35.588773 ], [ -72.203545, -35.590184 ], [ -72.200097, -35.59319  ], [ -72.194214, -35.586267 ], [ -72.193014, -35.581719 ], [ -72.18613,  -35.575239 ], [ -72.17897, -35.574595  ], [ -72.16969,  -35.584834 ], [ -72.173583, -35.593059 ], [ -72.176048, -35.60382  ], [ -72.18497,  -35.610247 ], [ -72.187715, -35.615171 ], [ -72.184774, -35.618164 ], [ -72.181208, -35.618258 ], [ -72.175342, -35.624659 ], [ -72.169229, -35.62482  ], [ -72.167881, -35.629436 ], [ -72.159844, -35.632562 ], [ -72.153286, -35.6344   ], [ -72.150229, -35.63448  ], [ -72.146939, -35.641645 ], [ -72.146132, -35.647079 ], [ -72.141742, -35.652191 ], [ -72.138321, -35.656028 ], [ -72.13812, -35.663945  ], [ -72.126176, -35.645517 ], [ -72.092,    -35.645568 ], [ -72.085535, -35.649898 ], [ -72.084659, -35.653668 ], [ -72.084836, -35.658244 ], [ -72.078497, -35.665902 ], [ -72.081174, -35.669165 ], [ -72.083835, -35.672011 ], [ -72.033842, -35.67245  ], [ -72.037521, -35.675272 ], [ -72.039224, -35.679809 ], [ -72.038951, -35.686062 ], [ -72.041117, -35.689339 ], [ -72.042773, -35.692628 ], [ -72.0429,   -35.695957 ], [ -72.047138, -35.700013 ], [ -72.050309, -35.702848 ], [ -72.055409, -35.702718 ], [ -72.058533, -35.704305 ], [ -72.061147, -35.705904 ], [ -72.059218, -35.708868 ], [ -72.065651, -35.711351 ],
    [ -72.067506, -35.712413 ], [ -72.071994, -35.711189 ], [ -72.073439, -35.713371 ], [ -72.077573, -35.714745 ], [ -72.079938, -35.717273 ], [ -72.082713, -35.718682 ], [ -72.087216, -35.717826 ], [ -72.091747, -35.71771  ], [ -72.095061, -35.721323 ], [ -72.099053, -35.719001 ], [ -72.101843, -35.720779 ], [ -72.103117, -35.718527 ], [ -72.108935, -35.716527 ], [ -72.112603, -35.717542 ], [ -72.114314, -35.714908 ], [ -72.118672, -35.710357 ], [ -72.122224, -35.708415 ], [ -72.128817, -35.703066 ], [ -72.133713, -35.700719 ], [ -72.138228, -35.700231 ], [ -72.141939, -35.702354 ], [ -72.143005, -35.706394 ], [ -72.146219, -35.70742  ], [ -72.145853, -35.709649 ], [ -72.150909, -35.711366 ], [ -72.155965, -35.713082 ], [ -72.161898, -35.714036 ], [ -72.166531, -35.716503 ], [ -72.175577, -35.715895 ], [ -72.183909, -35.720113 ], [ -72.191113, -35.718812 ], [ -72.197382, -35.716796 ], [ -72.201459, -35.716688 ], [ -72.204267, -35.718832 ], [ -72.202543, -35.721097 ], [ -72.204414, -35.722527 ], [ -72.222641, -35.724629 ], [ -72.227261, -35.726724 ], [ -72.235083, -35.729473 ], [ -72.238932, -35.734917 ], [ -72.24067,  -35.733021 ], [ -72.245489, -35.728822 ], [ -72.249357, -35.72354  ], [ -72.25428, -35.721927  ], [ -72.252846, -35.720117 ], [ -72.253193, -35.717518 ], [ -72.253028, -35.713454 ], [ -72.254297, -35.711201 ],
    [ -72.259189, -35.708849 ], [ -72.260352, -35.704009 ], [ -72.26704,  -35.701239 ], [ -72.269079, -35.695636 ], [ -72.272158, -35.693333 ], [ -72.276672, -35.69284  ], [ -72.281246, -35.693825 ], [ -72.286258, -35.694428 ], [ -72.289096, -35.697309 ], [ -72.296826, -35.697837 ], [ -72.299889, -35.695164 ], [ -72.303904, -35.693574 ], [ -72.30712,  -35.694595 ], [ -72.310712, -35.693757 ], [ -72.314619, -35.689581 ], [ -72.31495,  -35.686613 ], [ -72.317015, -35.681748 ], [ -72.32,     -35.677227 ], [ -72.322548, -35.673089 ], [ -72.32735, -35.674682  ], [ -72.33262,  -35.6787   ], [ -72.339437, -35.683091 ], [ -72.34615,  -35.684986 ], [ -72.349331, -35.687812 ],
    [ -72.355587, -35.690969 ], [ -72.362742, -35.691185 ], [ -72.366785, -35.69024  ], [ -72.371548, -35.69427  ], [ -72.375996, -35.690814 ], [ -72.383994, -35.686842 ], [ -72.391534, -35.684131 ], [ -72.394996, -35.681535 ], [ -72.398494, -35.679771 ], [ -72.403081, -35.679641 ], [ -72.409038, -35.675725 ], [ -72.411586, -35.675653 ], [ -72.418087, -35.672554 ], [ -72.421583, -35.670789 ], [ -72.426063, -35.668163 ], [ -72.428451, -35.664348 ], [ -72.435916, -35.659971 ], [ -72.437177, -35.65369  ], [ -72.438983, -35.648225 ], [ -72.444623, -35.648897 ], [ -72.45003,  -35.644162 ], [ -72.453614, -35.644475 ], [ -72.458199, -35.644344 ], [ -72.460111, -35.641374 ],
    [ -72.45944,  -35.637646 ], [ -72.454293, -35.636545 ], [ -72.450834, -35.639142 ], [ -72.447617, -35.635487 ], [ -72.444069, -35.636005 ], [ -72.44247,  -35.634386 ], [ -72.446366, -35.63011  ], [ -72.443712, -35.627688 ], [ -72.450503, -35.619583 ], [ -72.450413, -35.617503 ], [-72.437995,  -35.613279 ], [ -72.43729,  -35.608719 ], [ -72.432582, -35.605939 ], [-72.441554,  -35.57737  ], [ -72.437778, -35.572481 ], [ -72.437582, -35.567907 ], [ -72.436315, -35.562114 ], [ -72.437558, -35.555416 ], [ -72.437362, -35.550842 ], [ -72.435007, -35.543414 ], [ -72.435883, -35.540058 ]
  ]
]);

// USER: DEFINE TIME PERIODS
// Create multiple time periods to examine. Add as needed

// FUNCTION: To define time periods
var analysisPeriod = function(periodName, preStartDate, preEndDate, postStartDate, postEndDate){
  var period = ee.Dictionary({
    name: periodName,
    dates: ee.List([
      ee.List([ ee.String(preStartDate), ee.String(preEndDate)]),
      ee.List([ ee.String(postStartDate), ee.String(postEndDate)]) 
    ])
  })
  return(period)
}

// Creating time periods p1, p2, p3...
var p1 = analysisPeriod('p1', '2017-05-01', '2017-07-01', '2017-09-30', '2017-11-01')
var p2 = analysisPeriod('p2', '2018-05-01', '2018-07-01', '2018-09-30', '2018-11-01')
var p3 = analysisPeriod('p3', '2019-05-01', '2019-07-01', '2019-09-30', '2019-11-01')

// SELECT PERIODS FOR ANALYSIS: Comment out those unused or add as needed
var periods =  ee.List([
  p1,
  p2,
  p3
]);

// ----------------------------------------------------------------------------------
// MAPPING SETTINGS - Adding data to the map in the GEE code editor
// ----------------------------------------------------------------------------------

//  USER !!! Choose what you would like to see on the map!

// SETTINGS: ANALYSIS PERIODS TO MAP
// Strings indicate which analysis periods types will be added to the map. Comment out as desired
var periodsToMap = [
  1,
  2,
  3
]

// SETTINGS: IMAGE TYPES TO MAP
// Strings indicate which image types will be added to the map for each period. Comment out image types as desired
var imageTypesToMap = [
  'preFire_mosaic',
  'preFire_cm_mosaic',
  'postFire_mosaic',
  'postFire_cm_mosaic',
  'preNBR',
  'postNBR',
  'dNBR'
]

// ----------------------------------------------------------------------------------
// EXPORT SETTINGS 
// ----------------------------------------------------------------------------------

//  USER !!! Choose what data you would like to export (if anything)
// Folder where data will be stored (Google Drive)
var outputFolder = 'GEE'

// Periods to export (comment out as desired)
var periodsToExport = [
  1,
  2,
  3
] 
print('periods to export: ', periodsToExport)

// Image types to export (comment out as desired)
var imageTypesToExport = [
  'preFire_mosaic',
  'preFire_cm_mosaic',
  'postFire_mosaic',
  'postFire_cm_mosaic',
  'preNBR',
  'postNBR',
  'dNBR'
]
print('image types to export: ', imageTypesToExport)



// =======================================================================================================
// FUNCTIONS
// =======================================================================================================

// Cloud Mask Functions

// Mask clouds for Sentinel 2 imagery
var maskS2sr = function(image) {
    // Bits 10 and 11 are clouds and cirrus, respectively.
    var cloudBitMask = ee.Number(2).pow(10).int();
    var cirrusBitMask = ee.Number(2).pow(11).int();
    // Get the pixel QA band.
    var qa = image.select('QA60');
    // All flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
        .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
    // Return the masked image, scaled to TOA reflectance, without the QA bands.
    return image.updateMask(mask)
        .copyProperties(image, ["system:time_start"]);
}

// Function to mask clouds from the pixel quality band of Landsat 8 SR data.
var maskL8sr = function(image) {
    // Bits 3 and 5 are cloud shadow and cloud, respectively.
    var cloudShadowBitMask = 1 << 3;
    var cloudsBitMask = 1 << 5;
    var snowBitMask = 1 << 4;
    // Get the pixel QA band.
    var qa = image.select('pixel_qa');
    // All flags should be set to zero, indicating clear conditions.
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
        .and(qa.bitwiseAnd(snowBitMask).eq(0));
    // Return the masked image, scaled to TOA reflectance, without the QA bands.
    return image.updateMask(mask)
        .select("B[0-9]*")
        .copyProperties(image, ["system:time_start"]);
  }
  
  // dNBR function: Calculate dNBR as the difference between two NBR images
  // This function just calculates the difference between the pre- and post-fire images
  var calc_dNBR = function(preFire_dNBR_image, postFire_dNBR_image){
    // The result is called delta NBR or dNBR
    var dNBR_unscaled = preFire_dNBR_image.subtract(postFire_dNBR_image);
    // Scale product to USGS standards
    var dNBR_scaled = dNBR_unscaled.multiply(1000);
    return dNBR_scaled
  }

// ------------------------------------------------------------------------
// PERFORM dNBR FUNCTION

// This function calculates dNBR to identify burned areas (and a measure of burn severity)
// given a set of dates and an area of interest.

// Arguments:

// timePeriod - a nested list following this format: [[pre-fire start date, pre-fire end date],[post-fire start date, post-fire end date]]
// Dates must be strings in the format 'YYYY-MM-DD', and dates must follow a chronological order. 
// Dates should be selected depending on the period of interest during which burned areas are being identified
// For example, if looking for fires occurring in March 2018. The pre-fire image could be a mosaic of images
// from Jan & Feb 2018 and the post-fire image could be a mosaic of April & May 2018 images. Therefore,
// timePeriod = [['2018-01-01', '2018-02-28'], ['2018-04-01', '2018-05-31']]

// platform - 'S2' for Sentinel-2 or 'L8' for Landsat 8

// areaOfInterest - a user-defined polygon geometry, feature, or feature collection

// FUNCTION DEFINITION: burn severity with dNBR
var analysis_dNBR = function(timePeriod, platform, areaOfInterest){
  //*******************************************************
  //------------------ Processing user inputs ------------
  //print('platform:' , platform)
  //print('area of interest:' , aoi)
  
  // Date boundares for PRE-fire period
  // print('timePeriod input', timePeriod);
  var periodName = ee.String(timePeriod.get('name'))
  var prePeriod_bounds = ee.List(timePeriod.get('dates')).get(0);
  var prefire_start    = ee.List(prePeriod_bounds).get(0);       
  var prefire_end      = ee.List(prePeriod_bounds).get(1);
  // print('prePeriod bounds:', prePeriod_bounds,'prefire_start:',prefire_start,'prefire_end:', prefire_end);
  
  // Date boundaries for POST-fire period
  var postPeriod_bounds = ee.List(timePeriod.get('dates')).get(1);
  var postfire_start    = ee.List(postPeriod_bounds).get(0);       
  var postfire_end      = ee.List(postPeriod_bounds).get(1);
  // print('postPeriod bounds:', postPeriod_bounds,'postfire_start:',postfire_start,'postfire_end:', postfire_end);
  
  // Define the platform/satellite (Landsat = 'L8' | Sentinel2 = 'S2')
  // print('Platform:', platform);

  // Print satellite platform and dates to console
  if (platform == 'S2' | platform == 's2') {
      var ImCol = 'COPERNICUS/S2';
      var pl = 'Sentinel-2';
  } else {
      var ImCol = 'LANDSAT/LC08/C01/T1_SR';
      var pl = 'Landsat 8';
  }
  //print(ee.String('Satellite data selected for analysis: ').cat(pl));
  //print(ee.String('Fire incident period (fire occurred between): ').cat(prefire_end).cat(' and ').cat(postfire_start));

  // Location: Set study area as map center.
  var areaOfInterest = ee.FeatureCollection(areaOfInterest);
  Map.centerObject(areaOfInterest, 10);

  //*******************************************************
  // SELECT & PROCESS IMAGERY -----------------------------

  // LANDSAT 8 -------------------------------------
  
  // ACQUIRE IMAGE COLLECTION --------------
  
  var imagery = ee.ImageCollection(ImCol);

  // In the following lines imagery will be collected in an ImageCollection, depending on the
  // location of our study area, a given time frame and the ratio of cloud cover.
  var prefireImCol = ee.ImageCollection(imagery
      // Filter by dates.
      .filterDate(prefire_start, prefire_end)
      // Filter by location.
      .filterBounds(areaOfInterest));

  // Select all images that overlap with the study area from a given time frame 
  var postfireImCol = ee.ImageCollection(imagery
    // Filter by dates.
    .filterDate(postfire_start, postfire_end)
    // Filter by location.
    .filterBounds(areaOfInterest));
  
  // CHECK: Add the clipped images to the console on the right
  // print("Pre-fire Image Collection: ", prefireImCol); 
  // print("Post-fire Image Collection: ",
  //   ee. ImageCollection(postfireImCol).getInfo(0)
  // );
      
      
  // APPLY CLOUD AND SNOW MASK
  // Apply platform-specific cloud mask
  if (platform == 'S2' | platform == 's2') {
      var prefire_CM_ImCol = prefireImCol.map(maskS2sr)
      var postfire_CM_ImCol = postfireImCol.map(maskS2sr)
  } else {
      var prefire_CM_ImCol = prefireImCol.map(maskL8sr);
      var postfire_CM_ImCol = postfireImCol.map(maskL8sr);
  }

  //----------------------- Mosaic and clip images to study area -----------------------------
  // This is especially important, if the collections created above contain more than one image
  // (if it is only one, the mosaic() does not affect the imagery).   

  // Generate mosaic for the original imagery (not strictly necessary)
  var pre_mos = prefireImCol.mosaic().clip(areaOfInterest);
  var post_mos = postfireImCol.mosaic().clip(areaOfInterest);

  // Generate mosaic for the cloud- and snow-masked imagery
  var pre_cm_mos = prefire_CM_ImCol.mosaic().clip(areaOfInterest);
  var post_cm_mos = postfire_CM_ImCol.mosaic().clip(areaOfInterest);

  // CHECK:
  // Add the clipped images to the console on the right
  //print("Pre-fire True Color Image: ", pre_mos); 
  //print("Post-fire True Color Image: ", post_mos);

  //------------------ Calculate NBR for pre- and post-fire images ---------------------------

  // Apply platform-specific NBR = (NIR-SWIR2) / (NIR+SWIR2)
  if (platform == 'S2' | platform == 's2') {
      var preNBR =  pre_cm_mos.normalizedDifference (['B8', 'B12']);
      var postNBR = post_cm_mos.normalizedDifference(['B8', 'B12']);
  } else {
      var preNBR =  pre_cm_mos.normalizedDifference (['B5', 'B7']);
      var postNBR = post_cm_mos.normalizedDifference(['B5', 'B7']);
  }

  // FOR EXPORTING IMAGES FOR SENTINEL
  // .select(["B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11","QA10","QA20","QA60"]);

  
  // CHECK
  // Add the NBR images to the console on the right
  //print("Pre-fire Normalized Burn Ratio: ", preNBR); 
  //print("Post-fire Normalized Burn Ratio: ", postNBR);
  
  
  //------------------ Calculate difference between pre- and post-fire images ----------------

  var dNBR = calc_dNBR(preNBR, postNBR)
  
  // var NBR_multiband = dNBR.addBands(preNBR, postNBR)

  // CHECK: Add the difference image to the console on the right
  //print('Difference Normalized Burn Ratio: ', dNBR);
  
  //------------------ Generate output to return
  var outputImages = ee.ImageCollection([
    ee.Image(pre_mos)
      .set('system:id', periodName.cat(ee.String('_preFire_mosaic')))
      .set('system:time_start', ee.Date(prefire_start))
      .set('system:time_end',   ee.Date(prefire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(pre_cm_mos).toUint16()
      .set('system:id', periodName.cat(ee.String('_preFire_cm_mosaic')))
      .set('system:time_start', ee.Date(prefire_start))
      .set('system:time_end',   ee.Date(prefire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(post_mos)
      .set('system:id', periodName.cat(ee.String('_postFire_mosaic')))
      .set('system:time_start',ee.Date(postfire_start))
      .set('system:time_end',  ee.Date(postfire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(post_cm_mos).toUint16()
      .set('system:id', periodName.cat(ee.String('_postFire_cm_mosaic')))
      .set('system:time_start',ee.Date(postfire_start))
      .set('system:time_end',  ee.Date(postfire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(preNBR)
      .set('system:id', periodName.cat(ee.String('_preNBR')))
      .set('system:time_start', ee.Date(prefire_start))
      .set('system:time_end',   ee.Date(prefire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(postNBR)
      .set('system:id', periodName.cat(ee.String('_postNBR')))
      .set('system:time_start', ee.Date(postfire_start))
      .set('system:time_end',   ee.Date(postfire_end))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName),
    ee.Image(dNBR)
      .set('system:id', periodName.cat(ee.String('_dNBR')))
      .set('system:time_start', ee.Date(prefire_end))
      .set('system:time_end',   ee.Date(postfire_start))
      .set('firePeriod_start', ee.Date(prefire_end))
      .set('firePeriod_end', ee.Date(postfire_start))  
      .set('period', periodName)
  ]);
  
  // print(ee.String('Function complete for period: ').cat(periodName));
  return outputImages
}

// FUNCTION: Wrapper function to process multiple time periods by mapping over the list of periods
// Note: Only one AOI (containing 1+ polygons) and one platform (L8 or S2) can be processed at a time
var analysis_dNBR_multiple = function(timePeriod_input){
  return ee.ImageCollection((analysis_dNBR(ee.Dictionary(timePeriod_input), platform, aoi)));
};

// ==================================================================================
// RUN THE ANALYSIS 
// ==================================================================================

// GENERATE RESULTS 
// var result_onePeriod = analysis_dNBR(p1, 'S2', aoi)
// print('result - one period', result_onePeriod)
var results_list = periods.map(analysis_dNBR_multiple);
print('results:', results_list);
// print('results (type):', typeof results);

// MERGE COLLECTIONS:
// Merge different collections by iterating through them
var merge_col = function(nextCollection, previousCollection){
  return ee.ImageCollection(previousCollection).merge(ee.ImageCollection(nextCollection))
};
// empty collection to store images
var empty_col = ee.ImageCollection([])
var results_col = ee.ImageCollection(
      results_list.iterate(merge_col, empty_col)
  );
print('results collection:', results_col);


// ==================================================================================
// RESULTS
// ==================================================================================

// By period (choose one)
var p2_images = results_col
  .filterMetadata('period','equals', p2.get('name'))
print('p2 image collection:', p2_images)

// By image type (choose one, example given for dNBR)
var dNBR_images = results_col
  // .filterMetadata('system:id', 'contains', 'preFire_mosaic')
  // .filterMetadata('system:id', 'contains', 'postFire_mosaic')
  // .filterMetadata('system:id', 'contains', 'preFire_cm_mosaic')
  // .filterMetadata('system:id', 'contains', 'postFire_cm_mosaic')
  // .filterMetadata('system:id', 'contains', 'preNBR')
  // .filterMetadata('system:id', 'contains', 'postNBR')
  .filterMetadata('system:id', 'contains', 'dNBR')
print('dNBR image collection:', dNBR_images)

// By time period using reducers
// Example: This will produce an image of the max dNBR value for the months analyzed 2017-2019
var dNBR_2017_2019 = results_col
  .filterDate('2017-01-01', '2019-12-31')
  .filterMetadata('system:id', 'contains', 'dNBR')
  .reduce(ee.Reducer.max())
print('Max dNBR 2017-2019:', dNBR_2017_2019)

//  Ways to combine images with reducers
// var dNBR_max = dNBR_images.reduce(ee.Reducer.max())
// print('dNBR max', dNBR_max)

// var dNBR_mean = dNBR_images.reduce(ee.Reducer.mean())
// print('dNBR mean', dNBR_mean)

// var dNBR_min = dNBR_images.reduce(ee.Reducer.min() )
// print('dNBR min', dNBR_min)


// Make periods sequences for iteration
var nPeriods = periods.length()
print('n periods analyzed:', nPeriods)


//==========================================================================================
//                                    ADD LAYERS TO MAP
//==========================================================================================

//---------------------------------- True Color Imagery ------------------------------------

// Apply platform-specific visualization parameters for true color images
if (platform == 'S2' | platform == 's2') {
  var vis = {bands: ['B4', 'B3', 'B2'], max: 2000, gamma: 1.5};
} else {
  var vis = {bands: ['B4', 'B3', 'B2'], min: 0, max: 4000, gamma: 1.5};
}

//--------------------------- Greyscale -------------------------------

// Two-color palettes
var palette_grey = ['white', 'black'];
var palette_veg  = ['green', 'black'];
var palette_red  = ['white', 'red'];

//  Palette with the burn severity colors
var palette_burnSeverity =['7a8737', 'acbe4d', '0ae042', 'fff70b', 'ffaf38', 'ff641b', 'a41fd6', 'ffffff'];


//------------------------- Burn Ratio Product - Classification ----------------------------

// Define an SLD style of discrete intervals to apply to the image.
var sld_intervals =
  '<RasterSymbolizer>' +
    '<ColorMap type="intervals" extended="false" >' +
      '<ColorMapEntry color="#ffffff" quantity="-500" label="-500"/>' +
      '<ColorMapEntry color="#7a8737" quantity="-250" label="-250" />' +
      '<ColorMapEntry color="#acbe4d" quantity="-100" label="-100" />' +
      '<ColorMapEntry color="#0ae042" quantity="100" label="100" />' +
      '<ColorMapEntry color="#fff70b" quantity="270" label="270" />' +
      '<ColorMapEntry color="#ffaf38" quantity="440" label="440" />' +
      '<ColorMapEntry color="#ff641b" quantity="660" label="660" />' +
      '<ColorMapEntry color="#a41fd6" quantity="2000" label="2000" />' +
    '</ColorMap>' +
  '</RasterSymbolizer>';
  
// Seperate result into 8 burn severity classes
var thresholds = ee.Image([-1000, -251, -101, 99, 269, 439, 659, 2000]);

// Add the area of interest to the map
Map.addLayer(aoi,{},'AOI');
Map.centerObject(aoi, 11)

//---------------------------------- ADD LAYERS TO MAP -------------------------------------

// Using the settings defined at the beginning of the script, add images to the map
for (var i in periodsToMap) {
  var pdNumber = periodsToMap[i];
  var pdName = 'p' + pdNumber;
  var pdDates_list = periods.get(Number(i)).getInfo()['dates'];
  var pdFilename_end = platform + '_pre_' + pdDates_list[0][0] + '_' + 'post_' + pdDates_list[1][0] + '_' + pdDates_list[1][1];
  var pd_images = ee.ImageCollection(results_col)
    .filterMetadata('period','equals', pdName)
  for (var j in imageTypesToMap){
    var pd_image = ee.Image(ee.ImageCollection(pd_images).filterMetadata('system:id', 'contains', imageTypesToMap[j]).first());
    if (imageTypesToMap[j] == 'preFire_cm_mosaic' ){
      Map.addLayer(pd_image, vis, pdName + ': Pre-fire mosaic (cm)')
    }
    if (imageTypesToMap[j] == 'postFire_cm_mosaic') {
      Map.addLayer(pd_image, vis, pdName + ': Post-fire mosaic (cm)')
    }
    if (imageTypesToMap[j] == 'preNBR') {
      Map.addLayer(pd_image, {min: -1, max: 1, palette: palette_grey}, pdName + ': Pre-fire NBR (cm)')
    }
    if (imageTypesToMap[j] == 'postNBR') {
      Map.addLayer(pd_image, {min: -1, max: 1, palette: palette_grey}, pdName + ': Post-fire NBR (cm)')
    }
    if (imageTypesToMap[j] == 'dNBR') {
      // Map.addLayer( pd_image, {min: -500, max: 1000, palette: palette_red}, pdName + ': dNBR');
      Map.addLayer(pd_image.sldStyle(sld_intervals), {}, pdName + ': dNBR');
    }
  }
}


//==========================================================================================
//                                    EXPORT DATA
//==========================================================================================

// Set image scale based on platform (20 vs. 30 m)
var exportScale = {'S2': 20, 'L8': 30};

for (var i in periodsToExport) {
  var pdNumber = periodsToMap[i];
  var pdName = 'p' + pdNumber;
  var pdDates_list = periods.get(Number(i)).getInfo()['dates'];
  var pdFilename_end = platform + '_pre_' + pdDates_list[0][0] + '_' + pdDates_list[0][1] + '_post_' + pdDates_list[1][0] + '_' + pdDates_list[1][1];
  var pd_images = ee.ImageCollection(results_col).filterMetadata('period','equals', pdName)
  for (var j in imageTypesToExport){
    var pd_image = ee.Image(ee.ImageCollection(pd_images)
      .filterMetadata('system:id', 'contains', imageTypesToExport[j]).first());
    Export.image.toDrive({ image: pd_image,
      scale: exportScale[platform],
      region: aoi, maxPixels: 1e10,
      folder: outputFolder,
      fileNamePrefix: pdName + '_' + imageTypesToExport[j] + '_' + pdFilename_end,
      description: pdName + '_Export_' + imageTypesToExport[j]
    });
  }
}


