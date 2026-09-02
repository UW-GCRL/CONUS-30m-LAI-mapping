// ============================================================
// CONUS Monthly LAI (30 m) — Interactive Viewer, Time Series & Downloader
// ============================================================
// Data: projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m  (gap-filled, 276 monthly images)
//   Companions: ..._QA (gap-fill flags, 276)  |  ..._nvalid (valid obs months/yr, 23)
// Encoding: int16, scale 0.001 (LAI = value*0.001); NoData = -32768.
// Projection: NAD83 / Conus Albers (EPSG:5070), 30 m.
// Retrieval: Landsat 5/7/8/9 C2 L2 via the biome-stratified Random Forest of Kang et al. (2021);
//   MODIS LAI is the training target only. Gaps filled with a climatology-anchored temporal method.
// Code & docs: https://github.com/UW-GCRL/CONUS-30m-LAI-mapping
// ============================================================

var LAI    = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m');
var QA     = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m_QA');
var NVALID = ee.ImageCollection('projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m_nvalid');
var monthNames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];

// ---- Visualization presets per layer ----
var VIS = {
  LAI: {min: 0, max: 6,
        palette: ['#ffffcc','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#005a32'],
        labels: ['0','1','2','3','4','5','6'], name: 'LAI (m² m⁻²)'},
  QA:  {min: 0, max: 5,
        palette: ['#4d4d4d','#2c7bb6','#abd9e9','#fdae61','#d7191c','#bdbdbd'],
        labels: ['0 orig','1 interp','2 clim≥3y','3 clim1-2y','4 spatial','255 none'], name: 'QA code'},
  NVALID: {min: 0, max: 12,
        palette: ['#d7191c','#fdae61','#ffffbf','#a6d96a','#1a9641'],
        labels: ['0','3','6','9','12'], name: 'Valid months / yr'}
};

// ---- State ----
var state = {layer: 'LAI', year: 2020, month: 7, opacity: 1.0};
var currentLayer = null;
var marker = null;

// ============================================================
// UI widgets
// ============================================================
var title = ui.Label('CONUS Monthly LAI — 30 m, 2000–2022', {fontSize:'20px', fontWeight:'bold', margin:'10px 0 0 10px'});
var subtitle = ui.Label('Landsat-based LAI (Kang et al. 2021), gap-filled  |  int16 × 0.001 = LAI  |  EPSG:5070',
                        {fontSize:'12px', color:'#555', margin:'0 0 8px 10px'});

var layerSelect = ui.Select({items:['LAI','QA','NVALID'], value:'LAI', style:{margin:'2px 10px', stretch:'horizontal'}});
var years = []; for (var y=2000;y<=2022;y++) years.push(y.toString());
var yearSelect  = ui.Select({items:years, value:'2020', style:{margin:'2px 6px'}});
var monthSelect = ui.Select({items:monthNames, value:'Jul', style:{margin:'2px 6px'}});
var opacity = ui.Slider({min:0, max:1, value:1, step:0.05, style:{width:'150px', margin:'2px 10px'}});
var displayBtn = ui.Button({label:'🔄 Display', style:{margin:'4px 10px', stretch:'horizontal'}});

var info = ui.Label('', {fontSize:'12px', margin:'4px 10px', color:'#333'});
var hint = ui.Label('Tip: click anywhere on the map to plot the full LAI time series for that pixel.',
                    {fontSize:'11px', color:'#1a73e8', margin:'2px 10px'});

var legendTitle = ui.Label('', {fontWeight:'bold', margin:'8px 0 2px 10px'});
var legendPanel = ui.Panel({style:{margin:'0 10px 8px 10px'}});
var chartPanel  = ui.Panel({style:{margin:'4px 6px', minHeight:'220px'}});
chartPanel.add(ui.Label('(click the map for a time series)', {fontSize:'11px', color:'#999', margin:'6px 10px'}));

var assetLink = ui.Label('Open asset in Code Editor', {fontSize:'12px', color:'#1a73e8', margin:'4px 10px'});
assetLink.setUrl('https://code.earthengine.google.com/?asset=projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m');
var repoLink = ui.Label('Code & documentation on GitHub', {fontSize:'12px', color:'#1a73e8', margin:'2px 10px'});
repoLink.setUrl('https://github.com/UW-GCRL/CONUS-30m-LAI-mapping');
var cite = ui.Label('Cite: You, Kang & Chen (2026), Scientific Data. [DOI TBD]  |  Algorithm: Kang et al. (2021) RSE 258, 112383.',
                    {fontSize:'10px', color:'#888', margin:'6px 10px'});

// ---- Legend builder (dynamic per layer) ----
function buildLegend(kind) {
  var v = VIS[kind];
  legendTitle.setValue(v.name);
  var swatches = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  for (var i=0;i<v.palette.length;i++)
    swatches.add(ui.Label('', {backgroundColor:v.palette[i], padding:'8px 14px', border:'1px solid #ccc'}));
  var labs = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
  for (var j=0;j<v.labels.length;j++)
    labs.add(ui.Label(v.labels[j], {fontSize:'9px', padding:'0 4px', textAlign:'center'}));
  legendPanel.clear(); legendPanel.add(swatches); legendPanel.add(labs);
}

// ============================================================
// Data helpers
// ============================================================
function selectedImage() {
  var yr = parseInt(yearSelect.getValue(),10);
  var mo = monthNames.indexOf(monthSelect.getValue()) + 1;
  state.year = yr; state.month = mo; state.layer = layerSelect.getValue();
  if (state.layer === 'LAI') {
    return LAI.filter(ee.Filter.eq('year',yr)).filter(ee.Filter.eq('month',mo)).first()
             .multiply(0.001).rename('LAI').updateMask(ee.Image(0).eq(0));
  } else if (state.layer === 'QA') {
    return QA.filter(ee.Filter.eq('year',yr)).filter(ee.Filter.eq('month',mo)).first().rename('QA');
  } else { // NVALID (yearly)
    return NVALID.filter(ee.Filter.eq('year',yr)).first().rename('nvalid');
  }
}

function display() {
  var img = selectedImage();
  var kind = state.layer;
  if (kind === 'LAI') img = img.updateMask(img.gt(0)); // show valid vegetation only
  buildLegend(kind);
  if (currentLayer) Map.remove(currentLayer);
  var label = kind + ' ' + (kind==='NVALID' ? state.year
              : (monthNames[state.month-1] + ' ' + state.year));
  currentLayer = ui.Map.Layer(img, VIS[kind], label, true, state.opacity);
  Map.add(currentLayer);

  // scene stats
  var red = ee.Reducer.mean().combine(ee.Reducer.minMax(),'',true);
  img.reduceRegion({reducer:red, geometry:Map.getBounds(true), scale:1000, maxPixels:1e9, bestEffort:true})
     .evaluate(function(r){
        if(!r) return;
        var keys = Object.keys(r);
        var mk = keys.filter(function(k){return k.indexOf('_mean')>-1;})[0];
        var mn = keys.filter(function(k){return k.indexOf('_min')>-1;})[0];
        var mx = keys.filter(function(k){return k.indexOf('_max')>-1;})[0];
        var f = function(x){return (r[x]==null)?'N/A':r[x].toFixed(kind==='LAI'?2:0);};
        info.setValue('📊 '+label+' (in view)  |  mean '+f(mk)+'  min '+f(mn)+'  max '+f(mx));
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
    'Export.image.toDrive({\n'+
    '  image: ee.ImageCollection(\n'+
    '    "projects/ee-hyou34/assets/CONUS_Monthly_LAI_30m")\n'+
    '    .filter(ee.Filter.eq("year", '+yr+'))\n'+
    '    .filter(ee.Filter.eq("month", '+state.month+')).first(),\n'+
    '  description: "LAI_'+yr+'_'+mo+'",\n'+
    '  region: geometry, crs: "EPSG:5070", scale: 30, maxPixels: 1e13\n'+
    '});  // int16, LAI = value * 0.001');
}

// ============================================================
// Wire up + layout
// ============================================================
opacity.onSlide(function(v){ state.opacity=v; if(currentLayer) currentLayer.setOpacity(v); });
displayBtn.onClick(function(){ display(); updateDownload(); });
layerSelect.onChange(function(){
  // NVALID is yearly: month has no effect
  monthSelect.setDisabled(layerSelect.getValue()==='NVALID');
});

var controls = ui.Panel([
  ui.Label('Layer:', {fontWeight:'bold', margin:'4px 0 0 10px'}), layerSelect,
  ui.Panel([ui.Label('Year', {margin:'6px 0 0 10px'}), yearSelect,
            ui.Label('Month', {margin:'6px 0 0 6px'}), monthSelect],
           ui.Panel.Layout.flow('horizontal')),
  ui.Panel([ui.Label('Opacity', {margin:'6px 0 0 10px'}), opacity], ui.Panel.Layout.flow('horizontal')),
  displayBtn
]);

var panel = ui.Panel({
  widgets: [title, subtitle, controls, info, hint,
            legendTitle, legendPanel,
            ui.Label('Pixel time series', {fontWeight:'bold', margin:'8px 0 0 10px'}), chartPanel,
            dlTitle, dlCode, assetLink, repoLink, cite],
  style: {width:'400px', padding:'0'}
});
ui.root.insert(0, panel);

Map.setCenter(-96, 38, 4);
Map.setOptions('TERRAIN');
display();          // initial LAI layer
updateDownload();
