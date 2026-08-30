/***************************************
 * CONUS-30m-LAI: JUNE 2022 RE-EXPORT (single-month patch)
 *
 * Derived from lai_state_monthly_export.js (UW-GCRL/CONUS-30m-LAI-mapping).
 * Purpose: regenerate the missing/empty 2022 June CONUS LAI.
 *
 * CHANGES vs. original (all clearly marked with  // <<< JUNE2022):
 *   1. YEAR_TO_EXPORT = 2022
 *   2. MONTH_TO_EXPORT = 6  (export only June instead of all 12 months)
 *   3. RUN_ALL_STATES flag: enqueue all 49 CONUS states in one run
 *      (set false to keep the original one-state-per-run behaviour via STATE_INDEX)
 *   All science functions are UNCHANGED from the original.
 ***************************************/

/********** CONFIG **********/
var YEAR_TO_EXPORT    = 2022;                // <<< JUNE2022 (was 2000)
var MONTH_TO_EXPORT   = 6;                   // <<< JUNE2022 (June only)
var RUN_ALL_STATES    = true;                // <<< JUNE2022 (true = loop all CONUS states)
var NONVEG            = false;
var CLOUD_LIMIT       = 50;
var REDUCER           = 'median';
var EXPORT_TO         = 'Drive';
var ASSET_ROOT        = 'users/your_username/LAI_state_monthly';
var DRIVE_FOLDER      = 'CONUS_LAI';
var EXPORT_SCALE      = 30;
var EXPORT_CRS        = 'EPSG:5070';
var MAX_PIXELS        = 1e13;
var REGION_SIMPLIFY_M = 500;

// Used only when RUN_ALL_STATES = false (original one-state-per-run mode):
var STATE_INDEX         = 3;
var STATE_CODE_OVERRIDE = '';

var LAI_version = '0.2.0';

/********** SOURCES (States) **********/
var statesFC = ee.FeatureCollection('TIGER/2018/States');
var EXCLUDE = ee.List(['AK','HI','PR','VI','GU','MP','AS']);
var baseStates = statesFC.filter(ee.Filter.inList('STUSPS', EXCLUDE).not());

var stateCodes = ee.List(baseStates.aggregate_array('STUSPS')).sort().getInfo();
print('All CONUS state codes (sorted):', stateCodes);

// Current state; reassigned per-state when RUN_ALL_STATES = true.
var STATE_CODE = (STATE_CODE_OVERRIDE && STATE_CODE_OVERRIDE.length)
  ? STATE_CODE_OVERRIDE
  : stateCodes[Math.min(Math.max(STATE_INDEX, 0), stateCodes.length - 1)];

/********** LAI FUNCTIONS (UNCHANGED) **********/

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

var scaleLandsat = function(image) {
  return image.select(['green','red','nir','swir1'])
    .multiply(0.0000275).add(-0.2).divide(0.0001)
    .addBands(image.select('pixel_qa'));
};

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

var maskLandsat = function(image) {
  var qa = image.select('pixel_qa');
  var cloud  = qa.bitwiseAnd(1 << 3).neq(0);
  var shadow = qa.bitwiseAnd(1 << 4).neq(0);
  var water  = qa.bitwiseAnd(1 << 7).neq(0);
  return image.updateMask(cloud.not()).updateMask(shadow.not()).updateMask(water.not());
};

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

/********** DATA FETCH (UNCHANGED) **********/
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

/********** MONTHLY COMPOSITE & EXPORT (UNCHANGED logic) **********/

function validMaskFromQA(qa) {
  var bad_input = qa.bitwiseAnd(1).neq(0);
  var bad_lai   = qa.bitwiseAnd(2).neq(0);
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

/********** PER-STATE DRIVER (JUNE only)  <<< JUNE2022 **********/
var REGION_MAX_ERROR_M = 1000;

function runStateExports(stateCode) {
  STATE_CODE = stateCode;  // used inside buildMonthlyLAI's metadata

  var stateGeomProj = baseStates
    .filter(ee.Filter.eq('STUSPS', stateCode))
    .union(REGION_MAX_ERROR_M)
    .geometry(ee.ErrorMargin(REGION_MAX_ERROR_M))
    .buffer(0, REGION_MAX_ERROR_M)
    .transform(EXPORT_CRS, REGION_MAX_ERROR_M);

  var exportRegion = stateGeomProj
    .bounds(ee.ErrorMargin(REGION_MAX_ERROR_M), EXPORT_CRS)
    .transform(EXPORT_CRS, REGION_MAX_ERROR_M);

  var im = buildMonthlyLAI(YEAR_TO_EXPORT, MONTH_TO_EXPORT, exportRegion);   // <<< June only
  exportMonthlyImageForState(im, YEAR_TO_EXPORT, MONTH_TO_EXPORT, stateCode, exportRegion);
}

/********** RUN  <<< JUNE2022 **********/
print('Re-exporting', YEAR_TO_EXPORT, 'month', MONTH_TO_EXPORT,
      RUN_ALL_STATES ? '(ALL CONUS states)' : ('(state ' + STATE_CODE + ')'));

if (RUN_ALL_STATES) {
  // Enqueues one Export task per state (49). Start them from the Tasks tab.
  stateCodes.forEach(function(sc) { runStateExports(sc); });
} else {
  runStateExports(STATE_CODE);
}
