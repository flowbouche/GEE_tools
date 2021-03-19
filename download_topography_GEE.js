// GOOGLE EARTH ENGINE: Access and export topographic data

// Adam Bouché

// This script allows the user to input a polygon and export clipped topographic
// layers from the STRM digital elevation model (or other elevation date by changing the
// image stored as the dataset object)

// _______________________________________________________________________________________

// USER IMPORT GEOMETRY ASSETS, assign the name: 'aoi_input'

// SET UP AREA OF INTEREST (AOI)
// Set name of the area of interest
var areaName = 'testArea';

// Convert geom to feature
var aoi = ee.FeatureCollection(aoi_input);

// IMPORT DATA
// Import STRM data
var dataset = ee.Image('NASA/NASADEM_HGT/001');
// Select elevation
var dem = dataset.select('elevation')

// PROCESS DATA
// Clip strm to aoi
var dem = dem.clip(aoi);
print('Digital elevation model', dem);
// Create terrain products
var terrain = ee.Terrain.products(dem);
print('Terrain layer:', terrain);

  
// --------------------------------------------------------------------------------
// VISUALIZE ----
// Visualization params - Modify values for vis_dem to reflect the elevation range in the AOI
var vis_dem = {min: -100, max: 3000, palette: ['192F9C', 'D41000']};
var vis_slope = {min: 0, max: 90, palette: ['192F9C', 'D41000']};
var vis_aspect = {min: 0, max: 359, palette: ['192F9C', 'D41000']};
var vis_hs = {min: 50, max: 255, palette: ['192F9C', 'D41000']};

// Add to map
Map.addLayer(dem, vis_dem,'DEM', true);
Map.addLayer(terrain.select('elevation'), vis_dem,    'elev',       true);
Map.addLayer(terrain.select('slope'),     vis_slope,  'slope',      true);
Map.addLayer(terrain.select('aspect'),    vis_aspect, 'aspect',     true);
Map.addLayer(terrain.select('hillshade'), vis_hs,     'hillshade',  true);

Map.addLayer(aoi.draw({color: 'ffffff', strokeWidth: 3}), {},'AOI', true, 0.5);
Map.centerObject(aoi, 12);

// Export DEM
Export.image.toDrive({
  image: terrain.select('elevation'),
  folder: 'GEE',
  scale: 30,
  description: 'DEM_'.concat(areaName),
  region: aoi,
  maxPixels: 1e10
});

Export.image.toDrive({
  image: terrain.select('slope'),
  folder: 'GEE',
  scale: 30,
  description: 'slope_'.concat(areaName),
  region: aoi,
  maxPixels: 1e10
});

Export.image.toDrive({
  image: terrain.select("aspect"),
  folder: 'GEE',
  scale: 30,
  description: 'aspect_'.concat(areaName),
  region: aoi,
  maxPixels: 1e10
});

Export.image.toDrive({
  image: terrain.select("hillshade"),
  folder: 'GEE',
  scale: 30,
  description: 'hillshade_'.concat(areaName),
  region: aoi,
  maxPixels: 1e10
});