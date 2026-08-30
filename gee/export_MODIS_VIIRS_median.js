/*******************************************************************************
 * STAGE 1 (GEE): export monthly-MEDIAN MODIS & VIIRS LAI + a PFT map, at 500 m,
 * on the CONUS grid (EPSG:5070) — small rasters to download for LOCAL validation.
 *
 * Monthly MEDIAN (not mean) to match our product's monthly-median compositing.
 * Stored int16 x1000 (LAI x 1000, nodata -32768) — same encoding as our product,
 * so the local script reads all sources identically.
 *
 * Outputs (Drive/LAI_validation), one file per year, 12 monthly bands each:
 *   MODISmed_YYYY.tif   (MOD15A2H Terra, 2000-2022)
 *   VIIRSmed_YYYY.tif   (VNP15A2H, 2012-2022)
 *   PFT_500m.tif        (NLCD 2019 -> 8 biomes; 0 = non-veg/masked)
 * Restricted-mode: these are 500 m + median, so they are light; if needed, comment
 * out years you already exported.
 *******************************************************************************/
var MON = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
var YEARS = ee.List.sequence(2000, 2022).getInfo();
var DRIVE_FOLDER = 'LAI_validation';
var CRS = 'EPSG:5070', SCALE = 500;

var conus = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.inList('STUSPS', ['AK','HI','PR','VI','GU','MP','AS']).not())
  .geometry().dissolve(1000);

// monthly-MEDIAN LAI (MODLAND good quality), physical units -> int16 x1000
function monthlyMedian(collId, laiBand, year, mi) {
  var start = ee.Date.fromYMD(year, mi + 1, 1), end = start.advance(1, 'month');
  var col = ee.ImageCollection(collId).filterDate(start, end);
  var good = ee.ImageCollection(col.map(function(im) {
    var q = im.select('FparLai_QC').bitwiseAnd(1).eq(0);          // MODLAND_QC bit0 == 0
    var l = im.select(laiBand);
    return l.updateMask(l.lte(100)).updateMask(q).multiply(0.1);  // -> LAI
  }));
  var med = ee.Image(ee.Algorithms.If(col.size().gt(0), good.median(), ee.Image(0).selfMask()));
  return med.multiply(1000).round().toInt16().unmask(-32768).rename(MON[mi]);
}

function exportYear(collId, laiBand, prefix, year) {
  var bands = [];
  for (var mi = 0; mi < 12; mi++) bands.push(monthlyMedian(collId, laiBand, year, mi));
  var img = ee.Image.cat(bands).clip(conus);
  Export.image.toDrive({
    image: img, description: prefix + '_' + year, folder: DRIVE_FOLDER,
    region: conus, scale: SCALE, crs: CRS, maxPixels: 1e12,
    fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
  });
}

// MODIS (all years) and VIIRS (2012+)
YEARS.forEach(function(y) { exportYear('MODIS/061/MOD15A2H', 'Lai_500m', 'MODISmed', y); });
YEARS.filter(function(y){return y>=2012;}).forEach(function(y){ exportYear('NOAA/VIIRS/001/VNP15A2H','LAI','VIIRSmed', y); });

// YEAR-MATCHED PFT maps: one per NLCD epoch used in the original retrieval, so each
// image year is stratified by the biome the retrieval assumed (land cover changes over time).
// Epoch -> NLCD release (same as the generation nlcd_dict). Local script maps year -> epoch.
var PFT_EPOCHS = [[2001,'2019_REL'],[2004,'2019_REL'],[2006,'2019_REL'],[2008,'2019_REL'],
                  [2011,'2019_REL'],[2013,'2019_REL'],[2016,'2019_REL'],[2019,'2019_REL'],[2021,'2021_REL']];
function nlcdBiome(epochYear, release) {
  var lc = ee.ImageCollection('USGS/NLCD_RELEASES/' + release + '/NLCD')
    .filter(ee.Filter.calendarRange(epochYear, epochYear, 'year')).first().select('landcover');
  return lc.remap([41,42,43,52,71,81,82,90,95], [1,2,3,4,5,5,6,7,8]).unmask(0).toByte().rename('pft');
}
PFT_EPOCHS.forEach(function(e) {
  Export.image.toDrive({image: nlcdBiome(e[0], e[1]).clip(conus), description: 'PFT_' + e[0],
    folder: DRIVE_FOLDER, region: conus, scale: SCALE, crs: CRS, maxPixels: 1e12});
});

print('Exports queued: MODISmed_YYYY (23) + VIIRSmed_YYYY (11) + PFT_YYYY (9 epochs) = 43 tasks.',
      'Download all to K:\\Hangkai\\CONUS_LAI\\validation\\ then run validate_local.py.',
      'PFT: 1 decid,2 evergreen,3 mixed,4 shrub,5 grass/pasture,6 crop,7 woody-wet,8 herb-wet');
