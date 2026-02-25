/***************************************
 * CONUS-30m-LAI: Observation Density Analysis
 *
 * Computes the number of valid Landsat observations per pixel per month
 * across CONUS without running the LAI retrieval (much faster).
 * Used for Technical Validation — observation density statistics.
 *
 * Outputs:
 *   1. Map layers: obs_count for a selected month/year
 *   2. CONUS-wide seasonal summary chart (mean obs_count by month)
 *   3. Optional: export obs_count GeoTIFF for selected months
 ***************************************/

/********** CONFIG **********/
// All years in the dataset (2000–2022)
var YEARS          = [2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
                      2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,
                      2020,2021,2022];
var CLOUD_LIMIT    = 50;       // scene cloud threshold (%)
var EXPORT_MAPS    = false;    // set true to export obs_count GeoTIFFs
var EXPORT_SCALE   = 500;      // scale for export (use 500m to keep files small)
var EXPORT_CRS     = 'EPSG:5070';
var DRIVE_FOLDER   = 'CONUS_LAI_obscount';

/********** CONUS REGION **********/
var conus = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.inList('STUSPS',
    ee.List(['AK','HI','PR','VI','GU','MP','AS'])).not())
  .union(1000)
  .geometry(ee.ErrorMargin(1000));

/********** FUNCTIONS **********/

// QA mask: keep only clear land pixels (no cloud, shadow, water)
var maskLandsatQA = function(image) {
  var qa = image.select('QA_PIXEL');
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var water  = qa.bitwiseAnd(1 << 7).neq(0);
  return image.updateMask(cloud.not().and(shadow.not()).and(water.not()));
};

// Load all Landsat C02 L2 for a date range and region, apply QA mask
var getLandsatValid = function(start, end, region) {
  var makeValid = function(col) {
    return col
      .filterDate(start, end)
      .filterBounds(region)
      .filterMetadata('CLOUD_COVER', 'less_than', CLOUD_LIMIT)
      .map(maskLandsatQA)
      .select(['SR_B1'], ['valid']);   // keep any band as placeholder for counting
  };

  // L5 uses SR_B1; L7/L8/L9 also have SR_B1 (blue for L7/L8/L9)
  var L5 = makeValid(ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'));
  var L7 = makeValid(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'));
  var L8 = makeValid(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'));
  var L9 = makeValid(ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'));

  return L5.merge(L7).merge(L8).merge(L9);
};

// Count valid observations for a given year and month over a region
var getObsCount = function(year, month, region) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');
  var coll  = getLandsatValid(start.format('YYYY-MM-dd'), end.format('YYYY-MM-dd'), region);
  return coll.count()
    .rename('obs_count')
    .set('year', year)
    .set('month', month)
    .set('system:time_start', start.millis());
};

/********** 1. MAP VISUALIZATION **********/
// Show obs_count for a specific month/year
var SHOW_YEAR  = 2010;
var SHOW_MONTH = 7;   // July (peak growing season)

var obsJul = getObsCount(SHOW_YEAR, SHOW_MONTH, conus);
var obsJan = getObsCount(SHOW_YEAR, 1, conus);

var obsVis = {min: 0, max: 10, palette: ['#d73027','#f46d43','#fdae61','#fee090',
                                          '#ffffbf','#e0f3f8','#74add1','#313695']};
Map.centerObject(conus, 5);
Map.addLayer(obsJan.clip(conus), obsVis, 'obs_count Jan ' + SHOW_YEAR);
Map.addLayer(obsJul.clip(conus), obsVis, 'obs_count Jul ' + SHOW_YEAR, false);

/********** 2. SEASONAL SUMMARY CHART **********/
// For each month (1–12), compute mean obs_count across CONUS pixels.
// We sample at a coarse scale (5 km) to keep computation fast.

var SAMPLE_SCALE = 5000;  // 5 km sampling for chart statistics

var monthlyStats = ee.List.sequence(1, 12).map(function(m) {
  m = ee.Number(m);

  // Average across multiple years for a stable seasonal signal
  var annualCounts = ee.ImageCollection(
    ee.List(YEARS).map(function(y) {
      return getObsCount(ee.Number(y), m, conus);
    })
  );

  var meanCount = annualCounts.mean().rename('obs_count');

  // Reduce to CONUS mean and median at coarse scale
  var stats = meanCount.reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.median(), '', true)
             .combine(ee.Reducer.percentile([10, 25, 75, 90]), '', true),
    geometry: conus,
    scale: SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });

  return ee.Feature(null, stats.set('month', m));
});

var statsFC = ee.FeatureCollection(monthlyStats);

// Chart: mean obs_count by month
var chart = ui.Chart.feature.byFeature({
  features: statsFC,
  xProperty: 'month',
  yProperties: ['obs_count_mean', 'obs_count_median']
})
.setChartType('LineChart')
.setOptions({
  title: 'CONUS Mean Monthly Observation Count (Landsat, ' + YEARS.join('/') + ')',
  hAxis: {
    title: 'Month',
    ticks: [
      {v:1,f:'Jan'},{v:2,f:'Feb'},{v:3,f:'Mar'},{v:4,f:'Apr'},
      {v:5,f:'May'},{v:6,f:'Jun'},{v:7,f:'Jul'},{v:8,f:'Aug'},
      {v:9,f:'Sep'},{v:10,f:'Oct'},{v:11,f:'Nov'},{v:12,f:'Dec'}
    ]
  },
  vAxis: {title: 'Mean obs_count (valid pixels per month)'},
  series: {
    0: {color: '#e41a1c', lineWidth: 2, pointSize: 4, label: 'Mean'},
    1: {color: '#377eb8', lineWidth: 2, pointSize: 4, label: 'Median'}
  },
  legend: {position: 'top'}
});
print(chart);
print('Monthly stats table (export from chart for manuscript values):', statsFC);

/********** 3. REGIONAL BREAKDOWN **********/
// Compare obs_count in selected CONUS regions (for spatial context in manuscript)

var regionStats = function(regionName, regionGeom, month, year) {
  var obs = getObsCount(year, month, regionGeom);
  var stats = obs.reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.median(), '', true),
    geometry: regionGeom,
    scale: SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });
  return stats.set('region', regionName).set('month', month).set('year', year);
};

// Key US regions
var states = ee.FeatureCollection('TIGER/2018/States');
var pnwGeom  = states.filter(ee.Filter.inList('STUSPS', ['WA','OR'])).union(1000).geometry();
var seGeom   = states.filter(ee.Filter.inList('STUSPS', ['FL','GA','AL','MS','SC'])).union(1000).geometry();
var midwestGeom = states.filter(ee.Filter.inList('STUSPS', ['IA','IL','IN','OH','MO'])).union(1000).geometry();
var swGeom   = states.filter(ee.Filter.inList('STUSPS', ['AZ','NM','NV','UT'])).union(1000).geometry();

// Print regional stats for Jan and Jul (for manuscript text)
print('--- Regional obs_count (Jan ' + SHOW_YEAR + ') ---');
print('Pacific Northwest:', regionStats('PNW', pnwGeom, 1, SHOW_YEAR));
print('Southeast:',         regionStats('SE',  seGeom,  1, SHOW_YEAR));
print('Midwest:',           regionStats('MW',  midwestGeom, 1, SHOW_YEAR));
print('Southwest:',         regionStats('SW',  swGeom,  1, SHOW_YEAR));

print('--- Regional obs_count (Jul ' + SHOW_YEAR + ') ---');
print('Pacific Northwest:', regionStats('PNW', pnwGeom, 7, SHOW_YEAR));
print('Southeast:',         regionStats('SE',  seGeom,  7, SHOW_YEAR));
print('Midwest:',           regionStats('MW',  midwestGeom, 7, SHOW_YEAR));
print('Southwest:',         regionStats('SW',  swGeom,  7, SHOW_YEAR));

/********** 4. ZERO-OBSERVATION FRACTION **********/
// Fraction of CONUS pixels with obs_count == 0 per month (gap fraction)
// Useful for the manuscript statement about missing data

var gapStats = ee.List.sequence(1, 12).map(function(m) {
  m = ee.Number(m);
  var annualCounts = ee.ImageCollection(
    ee.List(YEARS).map(function(y) {
      return getObsCount(ee.Number(y), m, conus);
    })
  ).mean();

  // Fraction of pixels with < 1 valid observation on average
  var gapMask = annualCounts.lt(1);
  var gapFrac = gapMask.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: conus,
    scale: SAMPLE_SCALE,
    maxPixels: 1e9,
    bestEffort: true
  });
  return ee.Feature(null, gapFrac.set('month', m));
});

var gapChart = ui.Chart.feature.byFeature({
  features: ee.FeatureCollection(gapStats),
  xProperty: 'month',
  yProperties: ['obs_count']
})
.setChartType('ColumnChart')
.setOptions({
  title: 'Fraction of CONUS pixels with < 1 valid obs per month (avg over ' + YEARS.join('/') + ')',
  hAxis: {
    title: 'Month',
    ticks: [
      {v:1,f:'Jan'},{v:2,f:'Feb'},{v:3,f:'Mar'},{v:4,f:'Apr'},
      {v:5,f:'May'},{v:6,f:'Jun'},{v:7,f:'Jul'},{v:8,f:'Aug'},
      {v:9,f:'Sep'},{v:10,f:'Oct'},{v:11,f:'Nov'},{v:12,f:'Dec'}
    ]
  },
  vAxis: {title: 'Gap fraction (0–1)', minValue: 0, maxValue: 1},
  colors: ['#d73027']
});
print(gapChart);

/********** 5. OPTIONAL EXPORT **********/
// Export obs_count maps at 500m for all months of SHOW_YEAR
if (EXPORT_MAPS) {
  ee.List.sequence(1, 12).getInfo().forEach(function(m) {
    var mm  = (m < 10 ? '0' + m : '' + m);
    var obs = getObsCount(SHOW_YEAR, m, conus);
    Export.image.toDrive({
      image:           obs.toUint16(),
      description:     'obscount_CONUS_' + SHOW_YEAR + '_' + mm,
      fileNamePrefix:  'obscount_CONUS_' + SHOW_YEAR + '_' + mm,
      folder:          DRIVE_FOLDER,
      region:          conus.bounds(),
      scale:           EXPORT_SCALE,
      crs:             EXPORT_CRS,
      maxPixels:       1e13,
      fileFormat:      'GeoTIFF',
      formatOptions:   {cloudOptimized: true}
    });
  });
  print('Export tasks enqueued for year ' + SHOW_YEAR);
}
