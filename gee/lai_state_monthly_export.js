/***************************************
 * CONUS-30m-LAI: Monthly LAI Export (Google Earth Engine)
 *
 * Generates monthly Leaf Area Index (LAI) composites for a single US state
 * at 30-meter resolution using Landsat C02 L2 surface reflectance and a
 * biome-stratified Random Forest retrieval model (Kang et al. 2021, RSE).
 *
 * Reference:
 *   Kang, Y., et al. (2021). A data-driven approach to estimate leaf area
 *   index for Landsat images over the contiguous US. Remote Sensing of
 *   Environment, 258, 112383. https://doi.org/10.1016/j.rse.2021.112383
 *
 * Algorithm repository:
 *   https://github.com/yanghuikang/Landsat-LAI
 *
 * Usage:
 *   - Set YEAR_TO_EXPORT and STATE_INDEX (or STATE_CODE_OVERRIDE) per run.
 *   - Exports 12 monthly GeoTIFF files (one per month) to Google Drive
 *     or Earth Engine Assets.
 *   - Repeat for each state to cover all 48 contiguous US states + DC.
 *
 * Output bands:
 *   - LAI       : uint16, scale factor = 0.01, range 0–65535
 *   - obs_count : uint16, number of valid Landsat observations
 *
 * Output CRS  : EPSG:5070 (Albers Equal Area Conic, CONUS)
 * Output scale: 30 m
 ***************************************/

/********** CONFIG **********/
var YEAR_TO_EXPORT    = 2000;                // change per run
var NONVEG            = false;               // true to compute LAI on non-veg pixels too
var CLOUD_LIMIT       = 50;                  // scene cloud threshold (%)
var REDUCER           = 'median';            // 'median' or 'mean'
var EXPORT_TO         = 'Drive';             // 'Asset' or 'Drive'
var ASSET_ROOT        = 'users/your_username/LAI_state_monthly';
var DRIVE_FOLDER      = 'CONUS_LAI';
var EXPORT_SCALE      = 30;                  // meters
var EXPORT_CRS        = 'EPSG:5070';         // Albers Equal Area (CONUS)
var MAX_PIXELS        = 1e13;
var REGION_SIMPLIFY_M = 500;

// Process just ONE state per run:
// - Leave STATE_CODE_OVERRIDE = '' to use STATE_INDEX.
// - Otherwise set STATE_CODE_OVERRIDE = 'AZ' to run a specific state.
var STATE_INDEX         = 3;                // 0,1,2,... increment each run
var STATE_CODE_OVERRIDE = '';

var LAI_version = '0.2.0';

/********** SOURCES (States) **********/
var statesFC = ee.FeatureCollection('TIGER/2018/States');
var EXCLUDE = ee.List(['AK','HI','PR','VI','GU','MP','AS']);
var baseStates = statesFC.filter(ee.Filter.inList('STUSPS', EXCLUDE).not());

var stateCodes = ee.List(baseStates.aggregate_array('STUSPS')).sort().getInfo();
print('All CONUS state codes (sorted):', stateCodes);

var STATE_CODE = (STATE_CODE_OVERRIDE && STATE_CODE_OVERRIDE.length)
  ? STATE_CODE_OVERRIDE
  : stateCodes[Math.min(Math.max(STATE_INDEX, 0), stateCodes.length - 1)];
print('Selected state this run:', STATE_CODE);

/********** LAI FUNCTIONS **********/

var processLandsat = function(image) {
  var renamed_image  = renameLandsat(image);
  var rescaled_image = scaleLandsat(renamed_image);
  var masked_image   = maskLandsat(rescaled_image);
  var final_image    = getVIs(masked_image);

  var sunElevation = ee.Number(image.get('SUN_ELEVATION'));
  var sunAzimuth   = ee.Number(image.get('SUN_AZIMUTH'));
  var props = image.propertyNames();

  var solarZenith = ee.Number(ee.Algorithms.If(
    props.contains('SUN_ELEVATION'),
    ee.Number(90).subtract(sunElevation),
    ee.Algorithms.If(
      props.contains('SOLAR_ZENITH_ANGLE'),
      ee.Number(image.get('SOLAR_ZENITH_ANGLE')),
      45
    )
  ));
  var solarAzimuth = ee.Number(ee.Algorithms.If(
    props.contains('SUN_AZIMUTH'),
    sunAzimuth,
    ee.Algorithms.If(
      props.contains('SOLAR_AZIMUTH_ANGLE'),
      ee.Number(image.get('SOLAR_AZIMUTH_ANGLE')),
      180
    )
  ));

  final_image = final_image.copyProperties(image).set({
    'SOLAR_AZIMUTH_ANGLE': solarAzimuth,
    'SOLAR_ZENITH_ANGLE' : solarZenith
  });
  return ee.Image(final_image);
};

// Scale C02 SR to C01-like scale (0–10000)
var scaleLandsat = function(image) {
  return image.select(['green','red','nir','swir1'])
    .multiply(0.0000275).add(-0.2).divide(0.0001)
    .addBands(image.select('pixel_qa'));
};

// Rename Landsat bands (C02 L2); treat L9 same as L8
var renameLandsat = function(image) {
  var spacecraft = ee.String(image.get('SPACECRAFT_ID'));
  var spacecraft_no = ee.Number.parse(spacecraft.slice(8,9));
  var isL8or9 = spacecraft_no.gte(8);
  var from = ee.Algorithms.If(isL8or9,
    ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL'],
    ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL']);
  var to = ['blue','green','red','nir','swir1','swir2','pixel_qa'];
  return image.select(from, to).set('SPACECRAFT_NO', spacecraft_no);
};

// QA mask: remove cloud, shadow, water
var maskLandsat = function(image) {
  var qa = image.select('pixel_qa');
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var water  = qa.bitwiseAnd(1 << 7).neq(0);
  return image.updateMask(cloud.not()).updateMask(shadow.not()).updateMask(water.not());
};

// Compute NDVI and NDWI (features for RF model)
var getVIs = function(img) {
  var NDVI = img.expression('float((b("nir") - b("red"))) / (b("nir") + b("red"))');
  var NDWI = img.expression('float((b("nir") - b("swir1"))) / (b("nir") + b("swir1"))');
  return img.addBands(NDVI.select([0],['NDVI']))
            .addBands(NDWI.select([0],['NDWI']));
};

var getAffineTransform = function(image) {
  var projection = image.projection();
  var json = ee.Dictionary(ee.Algorithms.Describe(projection));
  return ee.List(json.get('transform'));
};

// LAI QA band (3 LSBits)
// bit0: input reflectance out of training range
// bit1: predicted LAI out of valid range (0–8)
// bit2: non-vegetation biome (biome=0)
var getLAIQA = function(landsat, sensor, lai) {
  var red_max=5100, green_max=5100, nir_max=7100, swir1_max=7100, lai_max=8;

  var data = ee.FeatureCollection('projects/ee-yanghuikang/assets/LAI/LAI_train_convex_hull_by_sensor_v0_1_1');
  var subset = data.filterMetadata('sensor','equals',sensor).sort('index');
  var hull_array = subset.aggregate_array('in_hull');
  var hull_array_reshape = ee.Array(hull_array).reshape([10,10,10,10]);

  var image_scaled = landsat.select('red').divide(red_max).multiply(10).floor().toInt()
    .addBands(landsat.select('green').divide(green_max).multiply(10).floor().toInt())
    .addBands(landsat.select('nir').divide(nir_max).multiply(10).floor().toInt())
    .addBands(landsat.select('swir1').divide(swir1_max).multiply(10).floor().toInt());

  var range_mask = landsat.select('red').gte(0)
    .and(landsat.select('red').lt(red_max))
    .and(landsat.select('green').gte(0))
    .and(landsat.select('green').lt(green_max))
    .and(landsat.select('nir').gte(0))
    .and(landsat.select('nir').lt(nir_max))
    .and(landsat.select('swir1').gte(0))
    .and(landsat.select('swir1').lt(swir1_max));

  var hull_image = image_scaled.select('red').multiply(0)
    .add(ee.Image(hull_array_reshape)).updateMask(range_mask);

  var in_mask = hull_image.arrayGet(
    image_scaled.select(['red','green','nir','swir1']).updateMask(range_mask)
  ).unmask(0).not().toByte();

  var out_mask = lai.gte(0).and(lai.lte(lai_max)).not().int();
  var biome_mask = landsat.select('biome2').eq(0).int();

  var qa_band = in_mask.bitwiseOr(out_mask.leftShift(1)).bitwiseOr(biome_mask.leftShift(2)).toByte();
  return qa_band.rename('QA');
};

// Map image year to nearest NLCD epoch
var getTrainImg = function(image) {
  var year = ee.Date(image.get('system:time_start')).get('year').format('%d');
  var nlcd_dict = {
    '1984':'2001','1985':'2001','1986':'2001','1987':'2001','1988':'2001','1989':'2001',
    '1990':'2001','1991':'2001','1992':'2001','1993':'2001','1994':'2001','1995':'2001',
    '1996':'2001','1997':'2001','1998':'2001','1999':'2001','2000':'2001','2001':'2001','2002':'2001',
    '2003':'2004','2004':'2004','2005':'2004',
    '2006':'2006','2007':'2006',
    '2008':'2008','2009':'2008',
    '2010':'2011','2011':'2011','2012':'2011',
    '2013':'2013','2014':'2013',
    '2015':'2016','2016':'2016','2017':'2016',
    '2018':'2019','2019':'2019','2020':'2019','2021':'2021','2022':'2021','2023':'2021','2024':'2021'
  };
  var nlcd_year = ee.Number.parse(ee.Dictionary(nlcd_dict).get(year));

  var addYearProp = function(img) {
    var y = ee.Date(img.get('system:time_start')).get('year');
    return img.set('year', y);
  };
  var nlcd_coll = ee.ImageCollection('USGS/NLCD_RELEASES/2019_REL/NLCD')
      .merge(ee.ImageCollection('USGS/NLCD_RELEASES/2021_REL/NLCD'))
      .map(addYearProp);
  var nlcd = nlcd_coll.filter(ee.Filter.eq('year', nlcd_year)).first();

  // Map NLCD classes to LAI biome types (0=non-veg, 1–8=vegetation biomes)
  var fromList = [11,12,21,22,23,24,31,41,42,43,51,52,71,72,73,74,81,82,90,95];
  var toList   = [0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 0, 4, 5, 0, 0, 0, 5, 6, 7, 8];
  var biome = ee.Image(nlcd).select('landcover').remap(fromList, toList).rename('biome2');

  image = processLandsat(image);

  var mask_img   = image.select(['pixel_qa'],['mask']).multiply(0);
  var sunZenith  = ee.Number(image.get('SOLAR_ZENITH_ANGLE'));
  var sunAzimuth = ee.Number(image.get('SOLAR_AZIMUTH_ANGLE'));

  image = image
    .addBands(biome)
    .addBands(mask_img.add(ee.Image.pixelLonLat()).select(['longitude'],['lon']))
    .addBands(mask_img.add(ee.Image.pixelLonLat()).select(['latitude'],['lat']))
    .addBands(ee.Image.constant(sunZenith).rename('sun_zenith').toFloat())
    .addBands(ee.Image.constant(sunAzimuth).rename('sun_azimuth').toFloat())
    .addBands(mask_img.add(1))
    .set('nlcd_year', nlcd_year);
  return image;
};

// Train biome-specific RF model
var getRFModel = function(sensor, biome) {
  var dir = 'projects/ee-yanghuikang/assets/LAI/LAI_train_sample_CONUS_full_v0_1_1';
  var train = ee.FeatureCollection(dir).filterMetadata('sensor','equals',sensor);
  if (biome > 0) train = train.filterMetadata('biome2','equals', biome);

  var features = ['red','green','nir','swir1','lat','lon','NDVI','NDWI','sun_zenith','sun_azimuth'];
  var rf = ee.Classifier.smileRandomForest({
    numberOfTrees: 100,
    minLeafPopulation: 50,
    variablesPerSplit: 5
  }).setOutputMode('REGRESSION').train({
    features: train, classProperty: 'MCD_LAI', inputProperties: features
  });
  return rf;
};

var getLAIforBiome = function(image, rf_model, biome) {
  return image.updateMask(image.select('biome2').eq(ee.Number(biome)))
              .classify(rf_model, 'LAI');
};

// Retrieve LAI for a single Landsat scene
var getLAIImage = function(image, nonveg) {
  nonveg = (nonveg === null || nonveg === undefined) ? false : nonveg;

  var satDigit = ee.String(image.get('SPACECRAFT_ID')).slice(8,9);
  var sensor = ee.Dictionary({'5':'LT05','7':'LE07','8':'LC08','9':'LC08'}).get(satDigit);

  var train_img = getTrainImg(image);

  var biomes = [1,2,3,4,5,6,7,8];
  if (nonveg) biomes = [0,1,2,3,4,5,6,7,8];

  var lai_img = train_img.select(['mask'],['LAI']).multiply(0).add(9999).double();
  biomes.forEach(function(b) {
    var biome = ee.Number(b);
    lai_img = lai_img.where(
      train_img.select('biome2').eq(biome),
      getLAIforBiome(train_img, getRFModel(ee.String(sensor), biome), biome)
    );
  });

  lai_img = lai_img.updateMask(lai_img.neq(9999));
  var qa = getLAIQA(train_img, ee.String(sensor), lai_img);

  lai_img = lai_img.rename('LAI').multiply(100).round().clamp(0,65535).uint16()
                   .addBands(qa.byte());

  return ee.Image(lai_img.copyProperties(image))
    .set('system:time_start', image.get('system:time_start'))
    .set('LAI_scale_factor', 0.01)
    .set('LAI_NLCD_year', train_img.get('nlcd_year'))
    .set('LAI_nonveg', nonveg)
    .set('LAI_version', LAI_version);
};

/********** DATA FETCH **********/
var getLandsat = function(start, end, region, cloud_limit) {
  cloud_limit = (cloud_limit == null) ? 50 : cloud_limit;

  var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .filterDate(start, end).filterBounds(region)
    .filterMetadata('CLOUD_COVER','less_than', cloud_limit);
  var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterDate(start, end).filterBounds(region)
    .filterMetadata('CLOUD_COVER','less_than', cloud_limit);
  var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterDate(start, end).filterBounds(region)
    .filterMetadata('CLOUD_COVER','less_than', cloud_limit);
  var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterDate(start, end).filterBounds(region)
    .filterMetadata('CLOUD_COVER','less_than', cloud_limit);

  return L5.merge(L7).merge(L8).merge(L9).sort('system:time_start');
};

/********** MONTHLY COMPOSITE & EXPORT **********/

function validMaskFromQA(qa) {
  var bad_input = qa.bitwiseAnd(1).neq(0);  // bit0: input out of range
  var bad_lai   = qa.bitwiseAnd(2).neq(0);  // bit1: LAI out of range
  return bad_input.or(bad_lai).not();
}

function buildMonthlyLAI(year, month, region) {
  var start = ee.Date.fromYMD(year, month, 1);
  var end   = start.advance(1, 'month');

  var lstack = getLandsat(start.format('YYYY-MM-dd'), end.format('YYYY-MM-dd'), region, CLOUD_LIMIT);

  var laiColl = lstack.map(function(img) {
    var laiImg = getLAIImage(img, NONVEG);
    var qa = laiImg.select('QA');
    var valid = validMaskFromQA(qa);
    return laiImg.select('LAI').updateMask(valid);
  });

  var collSize = laiColl.size();
  var empty = ee.Image(0).rename('LAI').updateMask(ee.Image(0));
  var safeColl = ee.ImageCollection(ee.Algorithms.If(
    collSize.gt(0), laiColl, ee.ImageCollection([empty])
  ));

  var laiMonthly = (REDUCER === 'mean') ? safeColl.mean() : safeColl.median();
  laiMonthly = laiMonthly.rename('LAI').toUint16();

  var count = safeColl.count().rename('obs_count').toUint16();

  var yyyymm = start.format('YYYYMM');
  var out = laiMonthly.addBands(count)
    .set({
      'system:index': yyyymm,
      'state': STATE_CODE,
      'year': year,
      'month': month,
      'start': start.millis(),
      'end': end.millis(),
      'LAI_scale_factor': 0.01,
      'composite_reducer': REDUCER,
      'nonveg': NONVEG,
      'cloud_limit': CLOUD_LIMIT,
      'LAI_version': LAI_version
    })
    .reproject({crs: EXPORT_CRS, scale: EXPORT_SCALE});

  return out;
}

function exportMonthlyImageForState(img, year, month, stateCode, region) {
  var mm   = (month < 10 ? '0' + month : '' + month);
  var name = 'LAI_' + stateCode + '_' + year + '_' + mm;

  if (EXPORT_TO === 'Asset') {
    Export.image.toAsset({
      image: img,
      description: name,
      assetId: ASSET_ROOT + '/' + name,
      region: region,
      scale: EXPORT_SCALE,
      crs: EXPORT_CRS,
      maxPixels: MAX_PIXELS
    });
  } else {
    Export.image.toDrive({
      image: img,
      description: name,
      fileNamePrefix: name,
      folder: DRIVE_FOLDER,
      region: region,
      scale: EXPORT_SCALE,
      crs: EXPORT_CRS,
      maxPixels: MAX_PIXELS,
      fileFormat: 'GeoTIFF',
      formatOptions: { cloudOptimized: true }
    });
  }
}

/********** DEFINE STATE REGION **********/
var REGION_MAX_ERROR_M = 1000;
var stateGeomProj = baseStates
  .filter(ee.Filter.eq('STUSPS', STATE_CODE))
  .union(REGION_MAX_ERROR_M)
  .geometry(ee.ErrorMargin(REGION_MAX_ERROR_M))
  .buffer(0, REGION_MAX_ERROR_M)
  .transform(EXPORT_CRS, REGION_MAX_ERROR_M);

var exportRegion = stateGeomProj
  .bounds(ee.ErrorMargin(REGION_MAX_ERROR_M), EXPORT_CRS)
  .transform(EXPORT_CRS, REGION_MAX_ERROR_M);

Map.centerObject(stateGeomProj, 6);
Map.addLayer(stateGeomProj, {color: 'red'}, 'State geometry (for reference)');
Map.addLayer(exportRegion,  {color: 'yellow'}, 'Export region (bounding rect)');
print('Exporting state:', STATE_CODE, 'year:', YEAR_TO_EXPORT);

/********** ENQUEUE 12 MONTHLY EXPORTS **********/
ee.List.sequence(1, 12).getInfo().forEach(function(m) {
  var im = buildMonthlyLAI(YEAR_TO_EXPORT, m, exportRegion);
  exportMonthlyImageForState(im, YEAR_TO_EXPORT, m, STATE_CODE, exportRegion);
});

// Sanity check: show January obs_count
var janImg = buildMonthlyLAI(YEAR_TO_EXPORT, 1, exportRegion);
Map.addLayer(janImg.select('obs_count'), {min:0, max:20}, 'obs_count (January)');
