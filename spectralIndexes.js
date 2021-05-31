// SPECTRAL INDIX CALCULATION

// Adam Bouché

// BAND DICTIONARIES - To easily select the right bands for raster algebra to calcule spectral indexes

var l8_bandDict = {
  "AEROSOL": "B1",
  "BLUE": "B2",
  "GREEN": "B3",
  "RED": "B4",
  "NIR": "B5",
  "SWIR1": "B6",
  "SWIR2": "B7",
  "PAN": "B8",
  "CIRRUS": "B9",
  "TIR1": "B10",
  "TIR2": "B11"
}

var s2_bandDict = {
  "BLUE": "B2",
  "GREEN": "B3",
  "RED": "B4",
  "REDE1": "B5",
  "REDE2": "B6",
  "REDE3": "B7",
  "NIR": "B8",
  "REDE4": "B8A",
  "WATERVAPOR": "B9",
  "SWIR1": "B11",
  "SWIR2": "B12",  
}

// Function used to select the band dictionary based on the sensor being used (Landsat 8 = 'L8', Sentinel 2 = 'S2')
var selectBandDictionary = function(sensor){
  if (sensor == 'S2' | sensor == 's2'){
    var bandDictToUse = s2_bandDict;
  };
  if (sensor == 'L8' | sensor == 'l8') { 
    var bandDictToUse = l8_bandDict;
  };
  return bandDictToUse
};

// Functions: Cloud/snow masks

// Function to mask clouds from the pixel quality band of Sentinel-2 SR data.
function maskS2sr(image) {
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
function maskL8sr(image) {
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

function cloudMask_imgCol(imgCol, sensorCode) {
  if (sensorCode == 'S2'){
    var imgCol_masked = imgCol.map(maskS2sr);
  }
  if (sensorCode == 'L8'){
    var imgCol_masked = imgCol.map(maskL8sr);
  }
  return imgCol_masked
};

// ---------------------------------------------------------------------------------------------------
// INDEXES

// NBRI (Normalized Burned Ratio Index)
var calcIndex_NBR = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float(NIR - SWIR2) / float ( NIR + SWIR2 )', {
    'NIR': imageIn.select(dictToUse['NIR']),
    'SWIR2': imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};


// EVI (Enhanced Vegetation Index)
var calcIndex_EVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)))', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE': imageIn.select(dictToUse['BLUE'])
  })
  return imageOut
};

// NDVI (Normalized Difference Vegetation Index)
var calcIndex_NDVI = function(imageIn, sensorCode) {
  var dictToUse = selectBandDictionary(sensorCode);
  var imageOut = imageIn.expression('float ((NIR - RED) / (NIR + RED))', {
    'NIR': imageIn.select(dictToUse['NIR']),
    'RED': imageIn.select(dictToUse['RED'])
  })
  return imageOut
};

// NDMI (Normalized Diff Moisture[water] Index)
var calcIndex_NDMI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('(NIR - SWIR1) / (NIR + SWIR1)', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'SWIR1':  imageIn.select(dictToUse['SWIR1'])
  })
  return imageOut
};

// GLI (Green Leaf Index)
var calcIndex_GLI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((GREEN - RED) + (GREEN - BLUE)) / ((2 * GREEN) + RED + BLUE))', {
    'GREEN':  imageIn.select(dictToUse['GREEN']),
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE'])
  })
  return imageOut
};

// SAVI (Soil Adjusted Vegetation Index)
var calcIndex_SAVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((NIR - RED) / (NIR + RED + L))* (1+L))', {
    'L':  0.5, // Mean vegetation cover (0-1)
    'NIR':  imageIn.select(dictToUse['NIR']),
    'RED':  imageIn.select(dictToUse['RED'])
  })
  return imageOut
};

// // GCI (Green Chlorophyll Index)
var calcIndex_GCI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((NIR) / (GREEN)) - 1)', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// // RGR (Red Green Ratio)
var calcIndex_RGR = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((RED) / (GREEN))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// SIPI (Structure Insensitive Pigment Index)
var calcIndex_SIPI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((NIR - BLUE) / (NIR - RED))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE']),
    'NIR':  imageIn.select(dictToUse['NIR'])
  })
  return imageOut
};

// ARVI (Atmospherically Resistant Vegetation Index)
var calcIndex_ARVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((NIR - (2 * RED) + BLUE) / (NIR + (2 * RED) + BLUE))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE']),
    'NIR':  imageIn.select(dictToUse['NIR'])
  })
  return imageOut
};

// NDTI (Normalized Difference Tillage Index) - Van Deventer et al., 1997)
var calcIndex_NDTI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((SWIR1 - SWIR2) / (SWIR1 + SWIR2))', {
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'SWIR2':  imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};

// CRC (Crop Residue Cover Index)
var calcIndex_CRC = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((SWIR1 - GREEN) / (SWIR1 + GREEN))', {
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// DFI (Dead Fuel Index) - Cao et al., 2010
var calcIndex_DFI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((1 - SWIR2 / SWIR1 ) * (RED / NIR))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'NIR':  imageIn.select(dictToUse['NIR']),
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'SWIR2':  imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};

// ---------------------------------------------------------------------------------------------------
// Example script use

// LOAD IMAGE COLLECTION
// Sensor - Choose S2 or L8
var sensor = 'S2';
// var sensor = 'L8';

// Choose area
var aoi = ohio;

// Set dates, max cloud cover
var maxCloudCover = 20;
var periodStart = '2020-03-15';
var periodEnd = '2020-05-30';

// GET COLLECTION

if (sensor == 'S2') {
  var imgCol_code = 'COPERNICUS/S2_SR'
}
if (sensor == 'L8') {
  var imgCol_code = 'LANDSAT/LC08/C01/T1_SR'
}
// Filter image collection
var imgCol = ee.ImageCollection(imgCol_code)
    .filterDate (ee.String(periodStart), ee.String(periodEnd))
    //.filterBounds (geometry)
    // .filterMetadata ('CLOUDY_PIXEL_PERCENTAGE', 'Less_Than', ee.Number(maxCloudCover));
// Mask clouds
var imgCol = cloudMask_imgCol(imgCol, sensor);
// Make mean image
var meanImage = imgCol
  // .mosaic()
  .mean()
  .clip(aoi);

// print(imgCol_mosaic)

// var meanImage = ee.Image(s2_collection.mean())
//   .mosaic()

// // NBRI (Normalized Burned Ratio Index)
var NBRI = calcIndex_NBR(meanImage, 'S2')

// NDVI (Normalized Difference Vegetation Index)
var NDVI = calcIndex_NDVI(meanImage, sensor);

// NDMI (Normalized Diff Moisture[water] Index)
var NDMI = calcIndex_NDMI(meanImage, sensor);

// EVI (Enhanced Vegetation Index)
var EVI = calcIndex_EVI(meanImage, sensor);

// GLI (Green Leaf Index)
var GLI = calcIndex_GLI(meanImage, sensor);

// SAVI (Soil Adjusted Vegetation Index)
var SAVI = calcIndex_SAVI(meanImage, sensor)

// GCI (Green Chlorophyll Index)
var GCI = calcIndex_GCI(meanImage, sensor);

// RGR (Red Green Ratio)
var RGR = calcIndex_RGR(meanImage, sensor);

// SIPI (Structure Insensitive Pigment Index)
var SIPI = calcIndex_SIPI(meanImage, sensor);

// ARVI (Atmospherically Resistant Vegetation Index)
var ARVI = calcIndex_ARVI (meanImage, sensor);

// NDTI (Normalized Difference Tillage Index)
var NDTI = calcIndex_NDTI(meanImage, 'S2')

// CRC (Crop Residue Cover Index)
var CRC = calcIndex_CRC(meanImage, 'S2')

// DFI (Dead Fuel Index)
var DFI = calcIndex_DFI(meanImage, 'S2')

// // Symbology (one for all images)
var symbology = {
  max: 1, min: 0,
  palette: ['#0000ff', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901', '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01', '012E01', '011D01', '011301']
};
      
// // Representacion de indices de meanImage
Map.addLayer (meanImage, {max: 4000.0, min: 0.0, gamma: 1.0, bands: ['B4','B3','B2']}, 'Image (RGB natural color)', 0);
Map.addLayer (NBRI, symbology, 'NBRI', 0);
Map.addLayer (NDVI, symbology, 'NDVI', 1);
Map.addLayer (NDMI, symbology, 'NDMI', 0);
Map.addLayer (GCI,  symbology, 'GCI',  0);
Map.addLayer (RGR,  symbology, 'RGR',  0);
Map.addLayer (GLI,  symbology, 'GLI',  0);
Map.addLayer (SAVI, symbology, 'SAVI', 0);
Map.addLayer (SIPI, symbology, 'SIPI', 0);
Map.addLayer (EVI,  symbology, 'EVI',  0);
Map.addLayer (ARVI, symbology, 'ARVI', 1);
Map.addLayer (NDTI, symbology, 'NDTI', 1);
Map.addLayer (CRC,  symbology, 'CRC',  1);
Map.addLayer (DFI,  symbology, 'DFI',  1);

Map.centerObject(aoi, 12); 