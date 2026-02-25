/***************************************
 * CONUS-30m-LAI: Observation Density — Export Results
 *
 * Exports all observation density statistics to Google Drive:
 *
 *   CSV outputs (small, fast):
 *     1. monthly_obs_stats.csv    — mean/median/p10/p25/p75/p90 by month
 *     2. regional_obs_stats.csv  — mean/median by region × month
 *     3. gap_fraction.csv        — % CONUS pixels with 0 valid obs by month
 *
 *   GeoTIFF outputs (one per month, 500 m resolution):
 *     4. obscount_CONUS_YYYY_MM.tif — mean obs_count maps
 *
 * How to use:
 *   1. Paste into GEE code editor and click Run.
 *   2. Go to the Tasks tab → click "Run" on each export task.
 *   3. Files appear in Google Drive under DRIVE_FOLDER.
 ***************************************/

/********** CONFIG **********/
var YEARS        = [2000, 2005, 2010, 2015, 2020, 2022]; // years to average over
var CLOUD_LIMIT  = 50;          // scene-level cloud cover threshold (%)
var SAMPLE_SCALE = 5000;        // 5 km — for CSV stats (fast)
var MAP_SCALE    = 500;         // 500 m — for exported GeoTIFF maps
var EXPORT_CRS   = 'EPSG:5070';
var DRIVE_FOLDER = 'CONUS_LAI_obscount';
var EXPORT_MAPS  = true;        // set false to skip GeoTIFF export

/********** CONUS GEOMETRY **********/
var conus = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.inList('STUSPS',
    ee.List(['AK','HI','PR','VI','GU','MP','AS'])).not())
  .union(1000)
  .geometry(ee.ErrorMargin(1000));

var conusBounds = conus.bounds();

/********** CORE FUNCTIONS **********/

var maskQA = function(image) {
  var qa     = image.select('QA_PIXEL');
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var water  = qa.bitwiseAnd(1 << 7).neq(0);
  return image.updateMask(cloud.not().and(shadow.not()).and(water.not()));
};

var getValidCollection = function(start, end, region) {
  var prep = function(col) {
    return col
      .filterDate(start, end)
      .filterBounds(region)
      .filterMetadata('CLOUD_COVER', 'less_than', CLOUD_LIMIT)
      .map(maskQA)
      .select(['SR_B1'], ['valid']);
  };
  return prep(ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'))
    .merge(prep(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')))
    .merge(prep(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')))
    .merge(prep(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')));
};

// Compute obs_count image for a single year+month
var getObsCount = function(year, month, region) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');
  return getValidCollection(
      start.format('YYYY-MM-dd'), end.format('YYYY-MM-dd'), region)
    .count()
    .rename('obs_count')
    .set({'year': year, 'month': month,
          'system:time_start': start.millis()});
};

// Average obs_count across multiple years for a given month
var getMeanObsCount = function(month, region) {
  var stack = ee.ImageCollection(
    ee.List(YEARS).map(function(y) {
      return getObsCount(ee.Number(y), ee.Number(month), region);
    })
  );
  return stack.mean().rename('obs_count');
};

/********** 1. MONTHLY STATS CSV **********/
// mean / median / p10 / p25 / p75 / p90 of obs_count across CONUS by month

var monthLabels = ['Jan','Feb','Mar','Apr','May','Jun',
                   'Jul','Aug','Sep','Oct','Nov','Dec'];

var monthlyStatsList = ee.List.sequence(1, 12).map(function(m) {
  var img   = getMeanObsCount(m, conus);
  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.median(),            '', true)
      .combine(ee.Reducer.percentile([10, 25, 75, 90]), '', true),
    geometry: conus,
    scale:    SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });
  var mm = ee.Number(m).int();
  return ee.Feature(null, stats
    .set('month_num',   mm)
    .set('month_label', ee.List(monthLabels).get(ee.Number(m).subtract(1)))
  );
});

var monthlyStatsFC = ee.FeatureCollection(monthlyStatsList);

Export.table.toDrive({
  collection:  monthlyStatsFC,
  description: 'monthly_obs_stats',
  folder:      DRIVE_FOLDER,
  fileNamePrefix: 'monthly_obs_stats',
  fileFormat:  'CSV',
  selectors:   ['month_num','month_label',
                'obs_count_mean','obs_count_median',
                'obs_count_p10','obs_count_p25',
                'obs_count_p75','obs_count_p90']
});
print('Export queued: monthly_obs_stats.csv');

/********** 2. REGIONAL STATS CSV **********/
// mean/median obs_count per region × month

var states = ee.FeatureCollection('TIGER/2018/States');

var regions = {
  'CONUS':      conus,
  'PacificNW':  states.filter(ee.Filter.inList('STUSPS', ['WA','OR'])).union(1000).geometry(),
  'Southeast':  states.filter(ee.Filter.inList('STUSPS', ['FL','GA','AL','MS','SC','NC'])).union(1000).geometry(),
  'Midwest':    states.filter(ee.Filter.inList('STUSPS', ['IA','IL','IN','OH','MO','MN','WI','MI'])).union(1000).geometry(),
  'Southwest':  states.filter(ee.Filter.inList('STUSPS', ['AZ','NM','NV','UT','CA'])).union(1000).geometry(),
  'Northeast':  states.filter(ee.Filter.inList('STUSPS', ['NY','PA','NJ','CT','MA','VT','NH','ME','RI'])).union(1000).geometry(),
  'GreatPlains':states.filter(ee.Filter.inList('STUSPS', ['KS','NE','SD','ND','OK','TX'])).union(1000).geometry()
};

var regionNames = Object.keys(regions);

var regionalStatsList = ee.List(regionNames).map(function(rName) {
  rName = ee.String(rName);
  var rGeom = ee.Dictionary(regions).get(rName);
  rGeom = ee.Geometry(rGeom);

  return ee.List.sequence(1, 12).map(function(m) {
    var img   = getMeanObsCount(m, rGeom);
    var stats = img.reduceRegion({
      reducer: ee.Reducer.mean()
        .combine(ee.Reducer.median(), '', true),
      geometry: rGeom,
      scale:    SAMPLE_SCALE,
      maxPixels: 1e9,
      bestEffort: true
    });
    return ee.Feature(null, stats
      .set('region',      rName)
      .set('month_num',   ee.Number(m).int())
      .set('month_label', ee.List(monthLabels).get(ee.Number(m).subtract(1)))
    );
  });
}).flatten();

var regionalStatsFC = ee.FeatureCollection(regionalStatsList);

Export.table.toDrive({
  collection:  regionalStatsFC,
  description: 'regional_obs_stats',
  folder:      DRIVE_FOLDER,
  fileNamePrefix: 'regional_obs_stats',
  fileFormat:  'CSV',
  selectors:   ['region','month_num','month_label',
                'obs_count_mean','obs_count_median']
});
print('Export queued: regional_obs_stats.csv');

/********** 3. GAP FRACTION CSV **********/
// Fraction of CONUS vegetated pixels with < 1 valid obs on average

// Load NLCD 2019 to mask out non-vegetated pixels for gap calculation
var nlcd = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')
  .filter(ee.Filter.calendarRange(2019, 2019, 'year')).first();
var vegMask = nlcd.select('landcover').remap(
  [11,12,21,22,23,24,31,41,42,43,51,52,71,72,73,74,81,82,90,95],
  [0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  0,  1,  1,  0,  0,  0,  1,  1,  1,  1]
).eq(1);  // 1 = vegetated pixel

var gapStatsList = ee.List.sequence(1, 12).map(function(m) {
  var img = getMeanObsCount(m, conus).updateMask(vegMask);

  // Fraction of veg pixels with < 1 obs on average
  var gapImg    = img.lt(1).rename('gap');
  var gapFrac   = gapImg.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: conus,
    scale:    SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });

  // Also compute fraction with < 2 obs
  var gap2Img   = img.lt(2).rename('gap_lt2');
  var gap2Frac  = gap2Img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: conus,
    scale:    SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });

  return ee.Feature(null, gapFrac
    .combine(gap2Frac)
    .set('month_num',   ee.Number(m).int())
    .set('month_label', ee.List(monthLabels).get(ee.Number(m).subtract(1)))
  );
});

var gapStatsFC = ee.FeatureCollection(gapStatsList);

Export.table.toDrive({
  collection:  gapStatsFC,
  description: 'gap_fraction',
  folder:      DRIVE_FOLDER,
  fileNamePrefix: 'gap_fraction',
  fileFormat:  'CSV',
  selectors:   ['month_num','month_label','gap','gap_lt2']
});
print('Export queued: gap_fraction.csv');

/********** 4. OBS_COUNT MAPS (GeoTIFF, 500 m) **********/
// Export mean obs_count maps (averaged across YEARS) for all 12 months

if (EXPORT_MAPS) {
  ee.List.sequence(1, 12).getInfo().forEach(function(m) {
    var mm  = (m < 10 ? '0' + m : '' + m);
    var img = getMeanObsCount(m, conus).toFloat();

    Export.image.toDrive({
      image:           img,
      description:     'obscount_map_' + mm,
      fileNamePrefix:  'obscount_CONUS_' + mm,
      folder:          DRIVE_FOLDER,
      region:          conusBounds,
      scale:           MAP_SCALE,
      crs:             EXPORT_CRS,
      maxPixels:       1e13,
      fileFormat:      'GeoTIFF',
      formatOptions:   {cloudOptimized: true}
    });
  });
  print('Export queued: 12 obs_count GeoTIFF maps (500 m).');
}

/********** QUICK MAP PREVIEW **********/
var obsVis = {
  min: 0, max: 8,
  palette: ['#d73027','#f46d43','#fdae61','#fee090',
            '#ffffbf','#74add1','#4575b4','#313695']
};
Map.centerObject(conus, 5);
Map.addLayer(getMeanObsCount(1, conus).clip(conus), obsVis, 'Mean obs_count Jan');
Map.addLayer(getMeanObsCount(7, conus).clip(conus), obsVis, 'Mean obs_count Jul', false);

print('=== All export tasks queued. Go to Tasks tab and click Run on each. ===');
print('Output folder in Google Drive:', DRIVE_FOLDER);
