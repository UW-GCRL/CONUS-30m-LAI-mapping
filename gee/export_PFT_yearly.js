/*******************************************************************************
 * Year-matched PFT (biome) maps for CONUS LAI validation, one file PER YEAR.
 *
 * MODIS LAI is already downloaded (MODISmed_YYYY.tif). This adds the missing PFT:
 * each year gets the NLCD epoch its LAI retrieval used, remapped to 8 biomes, at the
 * SAME 500 m EPSG:5070 grid as the MODIS files, so samples are stratified by the biome
 * that year's land cover (land cover changes across 2000-2022).
 *
 * Output (Drive/LAI_validation): PFT_2000.tif ... PFT_2022.tif  (23 files, 1 band, uint8)
 * PFT codes: 1 decid, 2 evergreen, 3 mixed, 4 shrub, 5 grass/pasture, 6 crop,
 *            7 woody wetland, 8 herbaceous wetland, 0 other/non-veg.
 * Download all 23 to:  K:\Hangkai\CONUS_LAI_gapfilled\MODIS\
 *******************************************************************************/
var DRIVE_FOLDER = 'LAI_validation';
var CRS = 'EPSG:5070', SCALE = 500;

var conus = ee.FeatureCollection('TIGER/2018/States')
  .filter(ee.Filter.inList('STUSPS', ['AK','HI','PR','VI','GU','MP','AS']).not())
  .geometry().dissolve(1000);

// year -> [NLCD epoch year, release]  (same nlcd_dict as the LAI generation)
var YEAR_EPOCH = {
  2000:[2001,'2019_REL'], 2001:[2001,'2019_REL'], 2002:[2001,'2019_REL'],
  2003:[2004,'2019_REL'], 2004:[2004,'2019_REL'], 2005:[2004,'2019_REL'],
  2006:[2006,'2019_REL'], 2007:[2006,'2019_REL'],
  2008:[2008,'2019_REL'], 2009:[2008,'2019_REL'],
  2010:[2011,'2019_REL'], 2011:[2011,'2019_REL'], 2012:[2011,'2019_REL'],
  2013:[2013,'2019_REL'], 2014:[2013,'2019_REL'],
  2015:[2016,'2019_REL'], 2016:[2016,'2019_REL'], 2017:[2016,'2019_REL'],
  2018:[2019,'2019_REL'], 2019:[2019,'2019_REL'], 2020:[2019,'2019_REL'],
  2021:[2021,'2021_REL'], 2022:[2021,'2021_REL']
};

function biome(epochYear, release) {
  var lc = ee.ImageCollection('USGS/NLCD_RELEASES/' + release + '/NLCD')
    .filter(ee.Filter.calendarRange(epochYear, epochYear, 'year')).first().select('landcover');
  //           41 42 43 52 71 81 82 90 95  ->  1 2 3 4 5 5 6 7 8
  return lc.remap([41,42,43,52,71,81,82,90,95], [1,2,3,4,5,5,6,7,8]).unmask(0).toByte().rename('pft');
}

Object.keys(YEAR_EPOCH).forEach(function(yStr) {
  var y = parseInt(yStr, 10), e = YEAR_EPOCH[y];
  Export.image.toDrive({
    image: biome(e[0], e[1]).clip(conus), description: 'PFT_' + y, folder: DRIVE_FOLDER,
    region: conus, scale: SCALE, crs: CRS, maxPixels: 1e12
  });
});

print('Queued 23 year-matched PFT exports: PFT_2000 ... PFT_2022.',
      'Run all tasks, download to K:\\Hangkai\\CONUS_LAI_gapfilled\\MODIS\\',
      'Codes: 1 decid,2 evergreen,3 mixed,4 shrub,5 grass/pasture,6 crop,7 woody-wet,8 herb-wet,0 other');
