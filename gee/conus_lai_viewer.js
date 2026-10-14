// ============================================================
// CONUS Monthly LAI (30 m) — Interactive Viewer, Time Series & Downloader
// ============================================================
// Data: projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m  (gap-filled, 276 monthly images)
//   Companions: ..._QA (gap-fill flags, 276)  |  ..._nvalid (valid obs months/yr, 23)
//   QA: 0 retrieved | 1 linear interpolation within the year | 2/3 climatology fill in a pixel-year
//   with >=3 / 1-2 valid months | 4 spatial-fallback climatology (<3 valid years for that month)
//   | 255 unfilled (no valid month in the pixel-year). Retrieved LAI = LAI where QA == 0.
// Encoding: int16, scale 0.001 (LAI = value*0.001); NoData = -32768 (masked; QA 255, nvalid 0).
// Grid: Albers Equal Area Conic, WGS 84 datum, EPSG:5070 parameters; 30 m, origin (-2417805, 3321225).
// Retrieval: Landsat 5/7/8/9 C2 L2 via the biome-stratified Random Forest of Kang et al. (2021);
//   MODIS LAI is the training target only. Gaps filled with a climatology-anchored temporal method.
// Viewer: the LAI, QA, and nvalid layers are drawn together (top to bottom); each has its own
//   on/off switch and opacity slider. Year and month apply to all layers (nvalid is yearly).
// Code & docs: https://github.com/UW-GCRL/CONUS-30m-LAI-mapping
// ============================================================

var LAI    = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m');
var QA     = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m_QA');
var NVALID = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m_nvalid');
var monthNames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];

// ---- Layers, drawn top to bottom in this order ----
var ORDER = ['LAI', 'QA', 'NVALID'];
var VIS = {
  LAI:    {min: 0, max: 6,
           palette: ['#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32']},
  QA:     {min: 0, max: 4,
           palette: ['#4d4d4d','#2c7bb6','#abd9e9','#fdae61','#d7191c']},
  NVALID: {min: 1, max: 12,
           palette: ['#d7191c','#fdae61','#ffffbf','#a6d96a','#1a9641']}
};
var LABEL = {LAI: 'LAI (m² m⁻²)', QA: 'Gap-fill QA', NVALID: 'nvalid (valid months / yr)'};
var cfg = {LAI:    {shown: true,  opacity: 1.0},
           QA:     {shown: false, opacity: 0.8},
           NVALID: {shown: false, opacity: 0.8}};
var mapLayers = {};

// ---- State ----
var state = {year: 2020, month: 7};
var marker = null;

// ============================================================
// UI widgets
// ============================================================
var title = ui.Label('CONUS Monthly LAI — 30 m, 2000–2022', {fontSize:'20px', fontWeight:'bold', margin:'10px 0 0 10px'});
var subtitle = ui.Label('Landsat-based LAI (Kang et al. 2021), gap-filled  |  int16 × 0.001 = LAI  |  30 m Albers (WGS 84)',
                        {fontSize:'12px', color:'#555', margin:'0 0 8px 10px'});

var years = []; for (var y=2000;y<=2022;y++) years.push(y.toString());
var yearSelect  = ui.Select({items:years, value:'2020', style:{margin:'2px 6px'}});
var monthSelect = ui.Select({items:monthNames, value:'Jul', style:{margin:'2px 6px'}});

var info = ui.Label('', {fontSize:'12px', margin:'4px 10px', color:'#333'});
var hint = ui.Label('Tip: click anywhere on the map to plot the full LAI time series for that pixel.',
                    {fontSize:'11px', color:'#1a73e8', margin:'2px 10px'});

var chartPanel  = ui.Panel({style:{margin:'4px 6px', minHeight:'220px'}});
chartPanel.add(ui.Label('(click the map for a time series)', {fontSize:'11px', color:'#999', margin:'6px 10px'}));

var assetLink = ui.Label('Open asset in Code Editor', {fontSize:'12px', color:'#1a73e8', margin:'4px 10px'});
assetLink.setUrl('https://code.earthengine.google.com/?asset=projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m');
var repoLink = ui.Label('Code & documentation on GitHub', {fontSize:'12px', color:'#1a73e8', margin:'2px 10px'});
repoLink.setUrl('https://github.com/UW-GCRL/CONUS-30m-LAI-mapping');
var cite = ui.Label('Data: You, Kang & Chen (2026), Zenodo, doi.org/10.5281/zenodo.23374236  |  Algorithm: Kang et al. (2021) RSE 258, 112383.',
                    {fontSize:'10px', color:'#888', margin:'6px 10px'});

// ---- Legends ----
function swatch(color) {
  return ui.Label('', {backgroundColor: color, width: '26px', padding: '6px 0', margin: '0',
                       border: '1px solid #ccc'});
}
function makeLegend(kind) {
  var p = ui.Panel({style: {margin: '0 10px 8px 34px'}});
  if (kind === 'QA') {
    var names = ['0 retrieved (not filled)', '1 interpolated within the year',
                 '2 climatology, year with ≥3 valid months', '3 climatology, year with 1–2 valid months',
                 '4 spatial-fallback climatology'];
    VIS.QA.palette.forEach(function(c, i) {
      p.add(ui.Panel([swatch(c), ui.Label(names[i], {fontSize: '11px', margin: '0 0 0 6px'})],
                     ui.Panel.Layout.flow('horizontal'), {margin: '1px 0'}));
    });
    p.add(ui.Label('255 (not filled) is transparent.', {fontSize: '10px', color: '#888', margin: '2px 0'}));
    return p;
  }
  var labels = (kind === 'LAI') ? ['0','1','2','3','4','5','6'] : ['1','','','','12'];
  var bar = ui.Panel([], ui.Panel.Layout.flow('horizontal'));
  var labs = ui.Panel([], ui.Panel.Layout.flow('horizontal'));
  VIS[kind].palette.forEach(function(c, i) {
    bar.add(swatch(c));
    labs.add(ui.Label(labels[i], {fontSize: '9px', width: '28px', textAlign: 'center', margin: '0'}));
  });
  p.add(bar); p.add(labs);
  return p;
}

// ---- Per-layer control: on/off + opacity + legend ----
function layerControl(kind) {
  var check = ui.Checkbox({label: LABEL[kind], value: cfg[kind].shown,
                           style: {fontWeight: 'bold', margin: '6px 0 0 10px'}});
  check.onChange(function(v) {
    cfg[kind].shown = v;
    if (mapLayers[kind]) mapLayers[kind].setShown(v);
  });
  var slider = ui.Slider({min: 0, max: 1, value: cfg[kind].opacity, step: 0.05,
                          style: {width: '110px', margin: '4px 10px 0 4px'}});
  slider.onSlide(function(v) {
    cfg[kind].opacity = v;
    if (mapLayers[kind]) mapLayers[kind].setOpacity(v);
  });
  var row = ui.Panel([check, ui.Label('opacity', {fontSize: '11px', color: '#666', margin: '9px 0 0 6px'}), slider],
                     ui.Panel.Layout.flow('horizontal'));
  return ui.Panel([row, makeLegend(kind)]);
}

// ============================================================
// Data helpers
// ============================================================
function imageFor(kind) {
  var yr = state.year, mo = state.month;
  if (kind === 'LAI') {
    var lai = LAI.filter(ee.Filter.eq('year', yr)).filter(ee.Filter.eq('month', mo)).first()
                 .multiply(0.001).rename('LAI');
    return lai.updateMask(lai.gt(0));
  }
  if (kind === 'QA') {
    return QA.filter(ee.Filter.eq('year', yr)).filter(ee.Filter.eq('month', mo)).first().rename('QA');
  }
  return NVALID.filter(ee.Filter.eq('year', yr)).first().rename('nvalid');   // yearly
}
function layerName(kind) {
  return (kind === 'NVALID') ? 'nvalid ' + state.year
                             : kind + ' ' + monthNames[state.month - 1] + ' ' + state.year;
}

// Redraw all layers for the selected year and month; on/off and opacity are kept.
function refresh() {
  state.year  = parseInt(yearSelect.getValue(), 10);
  state.month = monthNames.indexOf(monthSelect.getValue()) + 1;
  ORDER.slice().reverse().forEach(function(kind) {        // add bottom first, so LAI is on top
    var img = imageFor(kind);
    if (mapLayers[kind]) {
      mapLayers[kind].setEeObject(img);
      mapLayers[kind].setName(layerName(kind));
    } else {
      mapLayers[kind] = ui.Map.Layer(img, VIS[kind], layerName(kind), cfg[kind].shown, cfg[kind].opacity);
      Map.layers().add(mapLayers[kind]);
    }
  });
  updateStats();
  updateDownload();
}

// LAI statistics for the current view (computed on a 1 km sample)
function updateStats() {
  var name = layerName('LAI');
  info.setValue('📊 ' + name + ': computing...');
  imageFor('LAI').reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.minMax(), '', true),
    geometry: Map.getBounds(true), scale: 1000, maxPixels: 1e9, bestEffort: true
  }).evaluate(function(r) {
    if (!r || r.LAI_mean === null || r.LAI_mean === undefined) {
      info.setValue('📊 ' + name + ' (in view): no data'); return;
    }
    info.setValue('📊 ' + name + ' in view (1 km sample)  |  mean ' + r.LAI_mean.toFixed(2) +
                  '  min ' + r.LAI_min.toFixed(2) + '  max ' + r.LAI_max.toFixed(2));
  });
}

// ---- Click -> LAI time series at a pixel ----
function plotSeries(pt) {
  if (marker) Map.layers().remove(marker);
  marker = ui.Map.Layer(ee.FeatureCollection([ee.Feature(pt)]).style({color:'red', pointSize:6}), {}, 'clicked');
  Map.layers().add(marker);
  var col = LAI.map(function(img){
    var d = ee.Date.fromYMD(ee.Number(img.get('year')), ee.Number(img.get('month')), 1);
    return img.multiply(0.001).rename('LAI').set('system:time_start', d.millis());
  });
  var chart = ui.Chart.image.series({imageCollection: col, region: pt, reducer: ee.Reducer.mean(), scale: 30})
    .setChartType('LineChart')
    .setOptions({
      title: 'Monthly LAI, 2000–2022  ('+pt.coordinates().get(1).getInfo().toFixed(3)+'°, '+pt.coordinates().get(0).getInfo().toFixed(3)+'°)',
      hAxis:{title:'Date', format:'yyyy'}, vAxis:{title:'LAI (m² m⁻²)', viewWindow:{min:0}},
      lineWidth:1.2, pointSize:2, colors:['#238443'], legend:{position:'none'}, height:220
    });
  chartPanel.clear(); chartPanel.add(chart);
}
Map.onClick(function(coords){ plotSeries(ee.Geometry.Point([coords.lon, coords.lat])); });

// ---- Download code snippet ----
var dlTitle = ui.Label('Download (run in Code Editor):', {fontWeight:'bold', fontSize:'12px', margin:'8px 0 0 10px'});
var dlCode  = ui.Label('', {fontSize:'10px', margin:'2px 10px', whiteSpace:'pre'});
function updateDownload() {
  var yr = state.year, mo = monthNames[state.month-1];
  dlCode.setValue(
    'var img = ee.ImageCollection(\n'+
    '    "projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m")\n'+
    '    .filter(ee.Filter.eq("year", '+yr+'))\n'+
    '    .filter(ee.Filter.eq("month", '+state.month+')).first();\n'+
    'var grid = img.projection().getInfo();  // native 30 m grid\n'+
    'Export.image.toDrive({\n'+
    '  image: img.unmask(-32768), description: "LAI_'+yr+'_'+mo+'",\n'+
    '  region: geometry, crs: grid.wkt, crsTransform: grid.transform,\n'+
    '  maxPixels: 1e13, formatOptions: {noData: -32768}\n'+
    '});  // int16, LAI = value * 0.001; -32768 = no value');
}

// ============================================================
// Wire up + layout
// ============================================================
yearSelect.onChange(refresh);
monthSelect.onChange(refresh);

var dateRow = ui.Panel([ui.Label('Year', {margin:'8px 0 0 10px'}), yearSelect,
                        ui.Label('Month', {margin:'8px 0 0 6px'}), monthSelect],
                       ui.Panel.Layout.flow('horizontal'));
var layerPanel = ui.Panel([ui.Label('Layers (drawn top to bottom)', {fontWeight:'bold', margin:'8px 0 0 10px'})]
                          .concat(ORDER.map(layerControl)));

var panel = ui.Panel({
  widgets: [title, subtitle, dateRow, layerPanel, info, hint,
            ui.Label('Pixel time series', {fontWeight:'bold', margin:'8px 0 0 10px'}), chartPanel,
            dlTitle, dlCode, assetLink, repoLink, cite],
  style: {width:'400px', padding:'0'}
});
ui.root.insert(0, panel);

Map.setCenter(-96, 38, 4);
Map.setOptions('TERRAIN');
Map.setControlVisibility({layerList: false});   // layers are controlled from the side panel
refresh();          // draw all layers for the default year and month
