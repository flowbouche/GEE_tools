/*
Functions to perform cloud/snow masking and calculate spectral indexes
Currently designed for Landsat-8 and Sentinel-2 imagery

Author: Adam Bouché

To use these objects and functions, add the following code into your GEE script (change path and variable name as needed):
var spectral_indexes = require ('users/flow/modules:spectral_indexes')

Note: Mask functions compiled from GEE example scripts, spectral indexes compiled from scientific literature
*/


// MASK FUNCTIONS

// Function to mask clouds from the pixel quality band of Sentinel-2 SR data.
exports.mask_S2_sr = function(image) {
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
exports.mask_L8_sr = function(image) {
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

// MODIS: cloud and no data mask
// A function to mask out cloudy pixels and no data pixels.
exports.mask_modis = function(image) {
  // Cloud mask ----
  // Select the QA band.
  var QA = image.select('state_1km')
  // Make a mask to get bit 10, the internal_cloud_algorithm_flag bit.
  var bitMask = 1 << 10;
  // Make an image masking out cloudy areas.
  var image_couldMasked = image.updateMask(QA.bitwiseAnd(bitMask).eq(0));
  
  // No data mask ----
  // Find pixels that had observations.
  var withObs = image.select('num_observations_1km').gt(0);
  var image_masked = image_couldMasked.updateMask(withObs);
  
  return image_masked
}


// ---------------------------------------------------------------------------------------------------------



// BAND DICTIONARIES - To easily select the right bands for raster algebra to calcule spectral indexes
// Landsat-8 band dictionary
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

// Sentinel-2 band dictionary
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

// SPECTRAL INDEX FUNCTIONS

// NBRI (Normalized Burned Ratio Index)
exports.calcIndex_NBR = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float(NIR - SWIR2) / float ( NIR + SWIR2 )', {
    'NIR': imageIn.select(dictToUse['NIR']),
    'SWIR2': imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};


// EVI (Enhanced Vegetation Index)
exports.calcIndex_EVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)))', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE': imageIn.select(dictToUse['BLUE'])
  })
  return imageOut
};

// NDVI (Normalized Difference Vegetation Index)
exports.calcIndex_NDVI = function(imageIn, sensorCode) {
  var dictToUse = selectBandDictionary(sensorCode);
  var imageOut = imageIn.expression('float ((NIR - RED) / (NIR + RED))', {
    'NIR': imageIn.select(dictToUse['NIR']),
    'RED': imageIn.select(dictToUse['RED'])
  })
  return imageOut
};

// NDMI (Normalized Diff Moisture[water] Index)
exports.calcIndex_NDMI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('(NIR - SWIR1) / (NIR + SWIR1)', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'SWIR1':  imageIn.select(dictToUse['SWIR1'])
  })
  return imageOut
};

// GLI (Green Leaf Index)
exports.calcIndex_GLI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((GREEN - RED) + (GREEN - BLUE)) / ((2 * GREEN) + RED + BLUE))', {
    'GREEN':  imageIn.select(dictToUse['GREEN']),
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE'])
  })
  return imageOut
};

// SAVI (Soil Adjusted Vegetation Index)
exports.calcIndex_SAVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((NIR - RED) / (NIR + RED + L))* (1+L))', {
    'L':  0.5, // Mean vegetation cover (0-1)
    'NIR':  imageIn.select(dictToUse['NIR']),
    'RED':  imageIn.select(dictToUse['RED'])
  })
  return imageOut
};

// // GCI (Green Chlorophyll Index)
exports.calcIndex_GCI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float (((NIR) / (GREEN)) - 1)', {
    'NIR':  imageIn.select(dictToUse['NIR']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// // RGR (Red Green Ratio)
exports.calcIndex_RGR = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((RED) / (GREEN))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// SIPI (Structure Insensitive Pigment Index)
exports.calcIndex_SIPI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((NIR - BLUE) / (NIR - RED))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE']),
    'NIR':  imageIn.select(dictToUse['NIR'])
  })
  return imageOut
};

// ARVI (Atmospherically Resistant Vegetation Index)
exports.calcIndex_ARVI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((NIR - (2 * RED) + BLUE) / (NIR + (2 * RED) + BLUE))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'BLUE':  imageIn.select(dictToUse['BLUE']),
    'NIR':  imageIn.select(dictToUse['NIR'])
  })
  return imageOut
};

// NDTI (Normalized Difference Tillage Index) - Van Deventer et al., 1997)
exports.calcIndex_NDTI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((SWIR1 - SWIR2) / (SWIR1 + SWIR2))', {
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'SWIR2':  imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};

// CRC (Crop Residue Cover Index)
exports.calcIndex_CRC = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((SWIR1 - GREEN) / (SWIR1 + GREEN))', {
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'GREEN':  imageIn.select(dictToUse['GREEN'])
  })
  return imageOut
};

// DFI (Dead Fuel Index) - Cao et al., 2010
exports.calcIndex_DFI = function(imageIn, sensor) {
  var dictToUse = selectBandDictionary(sensor);
  var imageOut = imageIn.expression('float ((1 - SWIR2 / SWIR1 ) * (RED / NIR))', {
    'RED':  imageIn.select(dictToUse['RED']),
    'NIR':  imageIn.select(dictToUse['NIR']),
    'SWIR1':  imageIn.select(dictToUse['SWIR1']),
    'SWIR2':  imageIn.select(dictToUse['SWIR2'])
  })
  return imageOut
};
