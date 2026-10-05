/// Landsat Composite Generator

/// Imports

var Landsat5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2"),
    Landsat7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
    Landsat8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
    nirVis = {"bands":["SR_B5","SR_B4","SR_B3"],"min":0,"max":0.3},
    swir1vis = {"bands":["SR_B6","SR_B5","SR_B4"],"min":0,"max":0.3},
    Visible = {"bands":["SR_B4","SR_B3","SR_B2"],"min":0,"max":0.3},
    BoundingBox = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-79.21901142526966, 49.7711654622023],
          [-63.35574615642618, 49.7711654622023],
          [-63.35574615642618, 58.15441425790935],
          [-79.21901142526966, 58.15441425790935]]]),
    ShieldTaiga = ee.FeatureCollection("projects/ee-cesaramuradr1n24/assets/IntactForest_Quebec_ShieldTaiga"),
    provinceIFL = ee.FeatureCollection("projects/turnkey-topic-441614-c0/assets/EasternCanadianShieldTaigaIFL"),
    SampleArea = 
    /* color: #75ff16 */
    /* shown: false */
    ee.Geometry.Polygon(
        [[[-70.54159923775667, 55.98428411757292],
          [-70.54159923775667, 54.95368476234034],
          [-68.35531994088167, 54.95368476234034],
          [-68.35531994088167, 55.98428411757292]]], null, false),
    IFL_150 = ee.FeatureCollection("projects/ee-leafborealforests/assets/NAM_150_IFL"),
    NAM_IFL = ee.FeatureCollection("projects/ee-cesaramuradr1n24/assets/NAM_IFL"),
    GlobalWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater");

//////////////////////////////////// REQUIRED INPUTS ///////////////////////////////////////

var CLIP_EDGE = true;                 // Clip scene edges = true, false = keep whole scene
var EDGE_LENGTH_KM = 0.5;             // Length in kilometers to crop scenes edges
var BRIGHT_PIXEL_THRESHOLD  = 0.20;   // Mean reflectance threshold for saturated pixels
var CLOUD_COVER = 70;                 // Percentage (%) of cloud cover to filter scenes
var CALIBRATION = 'B';                // Q - 'Quebec Specific', B - 'Global Boreal' or N - 'None'
var COMPOSITE_METHOD = 'MaxNDVI';     // Options: 'MaxNDVI', 'Medoid', 'GeoMedian'
var LANDSAT_COLLECTION = 'A';         // Landsat collection A for L5+L7, B for L8+L7, C for L8 and D for L5
var STEP = 2;                         // Aggregation window
var START_YEAR = 2005;                // First composite year
var N_COMPOSITES = 1;                 // Number of composites
var POS_START = 6;                    // Peak of Season (POS) start month
var POS_END = 8;                      // Peak of Season (POS) end month
var WATER_MASK = false;               // Apply water mask: true = mask water pixels, false = keep them
var ROI = 'IFL_Q';                    // Select the region of interest (ROI): 'IFL_Q' (Quebec), 'IFL_C' (Custom) or 'Sample'
var IFL_ID = 'NAM_282_3';             // IFL ID corresponding to the IFL_C of interest

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// EXPORT CONTROL PANEL ///////////////////////////////////////

// Set to true/false to control which exports to run

var RUN_EXPORTS = {
  monthlyCounts:    false,    // total observations per month
  totalCounts:      false,    // total pixel observation count
  totalCoverage:    false,    // total spatial coverage
  landCoverage:     false,    // land-only spatial coverage
  compositeImages:  false,    // final composite GeoTIFFs
  meanNDVI:         false     // mean NDVI timeseries
};

var enabledExports = Object.keys(RUN_EXPORTS).filter(function(k) { return RUN_EXPORTS[k]; }).join('\n');

// State composite strategy for file export name convention
var compstrat = COMPOSITE_METHOD + '_' + CALIBRATION;

/////////////////////////////////////// COLLECTIONS TIME SPANS ////////////////////////////////////////

// Guard against mismatched collection and year choices
var validRanges = { 'A': [1984, 2011], 'B': [2013, 2024], 'C': [2013, 2024], 'D': [1984, 2011], 'E': [2013, 2024]};
var range = validRanges[LANDSAT_COLLECTION];

var END_YEAR = START_YEAR + (N_COMPOSITES * STEP) - 1;

if (START_YEAR < range[0] || END_YEAR > range[1]) {
  throw new Error(
    'Composite window ' + START_YEAR + '–' + END_YEAR + ' is outside the valid range for Collection ' 
    + LANDSAT_COLLECTION + ' (' + range[0] + '–' + range[1] + ')'
  );
}

////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// MASKING AND CALIBRAION //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

// Quality Assessment (QA) Masks for Landsat

function masksr(image) {
  var qa  = image.select('QA_PIXEL');
  var rad = image.select('QA_RADSAT');

  // Scale reflectance
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var output = image.addBands(opticalBands, null, true);

  // Basic mask
  var mask = qa.bitwiseAnd(1 << 0).eq(0)   // not fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))      // not dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))      // not cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))      // not cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // not shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0))      // not snow
    .and(rad.eq(0));                       // not saturated
  
  // Removes medium and high confidence cloud pixels
  // 0 = none, 1 = low, 2 = medium, 3 = high
  var cloudConf = qa.rightShift(8).bitwiseAnd(3);
  mask = mask.and(cloudConf.lte(1));        // keep only no/low confidence
  
  output = output.updateMask(mask);

  // Optional footprint edge clip
  if (CLIP_EDGE) {
    var smallGeom = ee.Geometry(output.geometry())
      .buffer(ee.Number(EDGE_LENGTH_KM).multiply(-1000));
    output = output.clip(smallGeom);
  }

  // require valid pixels across all reflectance bands
  var edgeMsk = output.select('SR_B.').mask().reduce(ee.Reducer.min());
  
  return output.updateMask(edgeMsk);
}

// Function to mask WATER pixels using GlobalWater dataset
function wtrmask(image) {
  var waterOccurrence = GlobalWater.select('occurrence');
  var fullWaterMask = waterOccurrence.gt(0).unmask(0).rename('FullWaterMask');
  var landMask = fullWaterMask.not();
  return image.updateMask(landMask);
}

// Conditional water mask — returns image unchanged if WATER_MASK is off
var applyWaterMask = function(image) {
  if (WATER_MASK) return wtrmask(image);
  return image;
};

var SR_BANDS = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];

var maskBrightPixels = function(img) {
  var bands = img.select(SR_BANDS);

  // Mean across all bands — elevated for saturated pixels
  var meanAllBands = bands.reduce(ee.Reducer.mean());

  var notSaturated = meanAllBands.lte(BRIGHT_PIXEL_THRESHOLD)
   
  return img.updateMask(notSaturated);
};

// Calculate NDVI from Landsat images
function computeNDVI(image) {
  var NDVI = image.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI').toFloat();
  return image.addBands(NDVI, null, true);
}


// Rename L5-L7 bands to match up L8
function rename(image){
  return image.select(
    ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7']);
}

//// Calibration Functions and Coefficients

// Derived calibration coefficients for Quebec
var L5coefficientsQ = {
  itcps: ee.Image.constant([-0.0006, 0.0008, -0.0028, 0.0185, 0.0069, 0.0023]), // Order: 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'
  slopes: ee.Image.constant([0.9048, 0.9000, 0.9637, 0.9149, 0.9452, 0.9560]),
  quadratic: ee.Image.constant([0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000]),
  cubic: ee.Image.constant([0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000])
};

var L8coefficientsQ = {
  itcps: ee.Image.constant([0.0109, 0.0075, 0.0049, 0.0276, 0.0077, 0.0061]), // Order: 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'
  slopes: ee.Image.constant([0.7828, 0.8714, 0.9777, 0.7384, 0.9353, 0.8393]),
  quadratic: ee.Image.constant([1.7087, 0.5558, 0.0000, 0.3850, 0.0000, 0.4892]),
  cubic: ee.Image.constant([0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000])
};

// // Berner et al. calibration coefficients
var L5coefficientsB = {
  itcps: ee.Image.constant([-0.0007, -0.0011, -0.0068, 0.0057, 0.0010, 0.0025]), // Order: 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'
  slopes: ee.Image.constant([0.8288, 0.8733, 1.0075, 0.9686, 0.9791, 0.9240]),
  quadratic: ee.Image.constant([1.8833, 0.3266, 0.0000, 0.0000, 0.0000, 0.3858]),
  cubic: ee.Image.constant([-6.0679, 0.0000, 0.0000, 0.0000, 0.0000, -0.5121])
};

var L8coefficientsB = {
  itcps: ee.Image.constant([0.0042, -0.0005, -0.0012, 0.0221, 0.0071, 0.0028]), // Order: 'Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2'
  slopes: ee.Image.constant([1.0950, 1.0412, 1.1592, 0.8442, 0.9698, 0.9212]),
  quadratic: ee.Image.constant([-0.4681, 0.0000, -0.9166, 0.1811, 0.0685, 0.5277]),
  cubic: ee.Image.constant([0.0000, 0.0000, 1.8795, 0.0000, 0.0000, -0.8244])
};

// Calibrate Landsat 5
function tmToETM(image) {
  var bands = image.select(SR_BANDS);
  var calibrateBands = bands.multiply(L5coefficients.slopes.toFloat())
    .add(L5coefficients.itcps.toFloat())
    .add(bands.pow(2).multiply(L5coefficients.quadratic.toFloat()))
    .add(bands.pow(3).multiply(L5coefficients.cubic.toFloat()))
    .toFloat();
  return image.addBands(calibrateBands, null, true)
}

// L7 is treated as the reference sensor, no calibration needed. This ensures structure matches.
function noCalibration(image) {
  return image.toFloat();
}

// Calibrate Landsat 8
function oliToETM(image) {
  var bands = image.select(SR_BANDS);
  var calibrateBands = bands.multiply(L8coefficients.slopes.toFloat())
    .add(L8coefficients.itcps.toFloat())
    .add(bands.pow(2).multiply(L8coefficients.quadratic.toFloat()))
    .add(bands.pow(3).multiply(L8coefficients.cubic.toFloat()))
    .toFloat();
  return image.addBands(calibrateBands, null, true)
}

// Calibration coefficient lookup: controlled by CALIBRATION switch
var calibrationOptions = {
  'Q':   { L5: L5coefficientsQ,  L8: L8coefficientsQ,  calFn: { L5: tmToETM,       L8: oliToETM      } },
  'B': { L5: L5coefficientsB, L8: L8coefficientsB,   calFn: { L5: tmToETM,       L8: oliToETM      } },
  'N':   { L5: null,             L8: null,             calFn: { L5: noCalibration, L8: noCalibration } }
};


var L5coefficients = calibrationOptions[CALIBRATION].L5;
var L8coefficients = calibrationOptions[CALIBRATION].L8;

var applyL5cal = calibrationOptions[CALIBRATION].calFn.L5;
var applyL8cal = calibrationOptions[CALIBRATION].calFn.L8;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// LOAD LANDSAT COLLECTIONS /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

var IFL_C = NAM_IFL.filter(ee.Filter.eq('IFL_ID', IFL_ID));

var roiOptions = {
  'IFL_Q': provinceIFL.geometry(),  // intact forest landscapes,
  'IFL_C': IFL_C.geometry(),
  'Sample': SampleArea              // small polygon for testing
};

var roi  = roiOptions[ROI];
var bbox = roi.bounds(); // Bounding Box

// Derived date range: do not edit these manually
var startYear = START_YEAR;
var endYear   = START_YEAR + (N_COMPOSITES * STEP) - 1;

// Date strings for filterDate (cover full years), summer filtering is handled by calendarRange separately
var startDate = ee.Date.fromYMD(startYear, 1, 1);
var endDate   = ee.Date.fromYMD(endYear,  12, 31);

// Reusable filter that applies study period and summer window
var filterCollection = function(collection) {
  return collection
    .filterBounds(bbox)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.calendarRange(POS_START, POS_END, 'month'));
};

// Fixed Landsat 7 SLC-off failure date
var SLC_FAILURE = '2003-06-01';

var landsat5 = filterCollection(Landsat5.filterDate('1984-03-16', '2011-12-31')) // Total range "1984-03-16","2012-05-05"
                      .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_COVER))
                      .map(masksr)  // Apply cloud and shadow mask
                      .map(applyWaterMask) // Apply water mask
                      .map(rename) // Must precede calibration
                      .map(applyL5cal);

var landsat7_slcon = filterCollection(Landsat7.filterDate('1999-05-28', SLC_FAILURE)) // Total range '1999-01-01','2003-05-30'
                        .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_COVER))
                        .map(masksr)  // Apply cloud and shadow mask
                        .map(applyWaterMask) // Apply water mask
                        .map(rename)
                        .map(noCalibration);
                        
var landsat7_slcoff_A = filterCollection(Landsat7.filterDate(SLC_FAILURE, '2011-12-31')) // Total range '2003-06-01','2024-01-19'
                        .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_COVER))
                        .map(masksr)  // Apply cloud and shadow mask
                        .map(applyWaterMask) // Apply water mask
                        .map(rename)
                        .map(noCalibration);

var landsat7_slcoff_B = filterCollection(Landsat7.filterDate('2013-01-01', '2024-01-19')) // Total range '2003-06-01','2024-01-19'
                        .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_COVER))
                        .map(masksr)  // Apply cloud and shadow mask
                        .map(applyWaterMask) // Apply water mask
                        .map(rename)
                        .map(noCalibration);

var landsat8 = filterCollection(Landsat8.filterDate('2013-03-18', '2024-12-31')) // Total range '2013-03-18' - present
                      .filter(ee.Filter.lte('CLOUD_COVER', CLOUD_COVER))
                      .map(masksr)  // Apply cloud and shadow mask
                      .map(applyWaterMask) // Apply water mask
                      .select(SR_BANDS)
                      .map(applyL8cal);

/////// Merge Landsat collections ////////
var collectionOptions = {
  // Landsat 5 + L7 SLC-on + L7 SLC-off: 1984-2011
  'A': landsat5.merge(landsat7_slcon).merge(landsat7_slcoff_A),

  // Landsat 8 + L7 SLC-off: 2013-present
  'B': landsat8.merge(landsat7_slcoff_B),
  
  // Landsat 8 ONLY: 2013-present
  'C': landsat8,
  
  // Landsat 5 ONLY: 1984-2011
  'D': landsat5
};

var landsat_collection_merged = collectionOptions[LANDSAT_COLLECTION];

// Function to get unique sensor names
function getSensors(collection) {
  return ee.List(collection.aggregate_histogram('SPACECRAFT_ID').keys());
}

var sensors = getSensors(landsat_collection_merged).join('_');
var expectedSensors = { 'A': ['LANDSAT_5', 'LANDSAT_7'], 'B': ['LANDSAT_7', 'LANDSAT_8'], 
  'C': ['LANDSAT_8'] , 'D': ['LANDSAT_5']};
var sensorWindows = {                                     // Approximate operational windows
  'LANDSAT_5': {start: '1984-03-16', end: '2012-05-05'},
  'LANDSAT_7': {start: '1999-05-28', end: '2024-01-19'},  
  'LANDSAT_8': {start: '2013-03-18', end: '2030-01-01'}
};
var rangeStart = new Date(startDate).getTime();
var rangeEnd   = new Date(endDate).getTime();
var expected = expectedSensors[LANDSAT_COLLECTION].filter(function(sensor) {
  var w = sensorWindows[sensor];
  var sensorStart = new Date(w.start).getTime();
  var sensorEnd   = new Date(w.end).getTime();
  return rangeStart <= sensorEnd && rangeEnd >= sensorStart;
});

sensors.evaluate(function(s) {
  expected.forEach(function(sensor) {
    if (s.indexOf(sensor) === -1) {
      throw new Error('Expected sensor ' + sensor + ' not found in Collection ' + LANDSAT_COLLECTION + 
        '. Check your date range and collection settings.');
    }
  });
});

var firstImage = landsat_collection_merged.sort('system:time_start').first(); 
var firstDate = ee.Date(firstImage.get('system:time_start')); // Get the earliest image (start date)
var firstYear = firstDate.get('year');
var firstMonth = firstDate.format('MM').getInfo();

var lastImage = landsat_collection_merged.sort('system:time_start', false).first();
var lastDate = ee.Date(lastImage.get('system:time_start')); // Get the latest image (end date)
var lastYear = lastDate.get('year');
var lastMonth = lastDate.format('MM').getInfo();

// Print parameters in console to check

print('SELECTED PARAMETERS', 'Calibration ' + CALIBRATION + ' using ' + COMPOSITE_METHOD + ' | Water mask: ' + (WATER_MASK ? 'ON' : 'OFF') + ' | ROI: ' + ROI);
print('COMPOSITE DETAILS', STEP + '-year composites' + ' from ' + sensors.getInfo())
print('Start Date:', firstDate.format('YYYY-MM-dd'));
print('End Date:', lastDate.format('YYYY-MM-dd'));
print('Exports enabled: ', (enabledExports || 'none'));

///////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Multi-year COMPOSITES generation /////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

// Define composite strategy functions

// Maximum NDVI composite: selects peak greenness pixel per compositing window
var maxNDVIComposite = function(collection) {
  return collection.map(computeNDVI)          // temporary NDVI for selection only
    .qualityMosaic('NDVI')
    .select(SR_BANDS);                       // drop NDVI, added cleanly later;
};

// Geometric Median composite: synthetic pixel value representative of the multivariate median across all bands
var geoMedianComposite = function(collection) {
  var nBands = SR_BANDS.length; // Explicitly select SR_BANDS, excludes any other bands
  return collection
    .select(SR_BANDS)
    .map(maskBrightPixels)
    .reduce(ee.Reducer.geometricMedian(nBands)) // Requires the number of bands explicitly
    .rename(SR_BANDS)  // Output bands are renamed with '_geometricMedian' suffix, rename back
    .toFloat();
};

// Change COMPOSITE_METHOD switch to swap approach
var makeComposite = function(collection) {
  if (COMPOSITE_METHOD === 'GeoMedian')        return geoMedianComposite(collection);
  else                                         return maxNDVIComposite(collection);
};

//// Build composites from merged collection
var multiYearComposite = function(startYear) {
  var start      = ee.Date.fromYMD(startYear, 1, 1);
  var end        = start.advance(STEP, 'year').advance(-1, 'day');
  var collection = landsat_collection_merged.filterDate(start, end)
  var composite = makeComposite(collection);

  return composite
    .addBands(computeNDVI(composite).select('NDVI'))
    .set('start_year', startYear)
    .set('end_year',   end.get('year'))
    .clip(roi)
};

var multiYearIntervals = ee.List.sequence(startYear, endYear, STEP);                      // Step by # years
var multiYearCollection = ee.ImageCollection(multiYearIntervals.map(multiYearComposite)); // Create a composite collection

var labels = [];
for (var i = 0; i < N_COMPOSITES; i++) {
  var compStart = START_YEAR + (i * STEP);
  var compEnd   = compStart + STEP - 1;
  labels.push('Composite ' + (i + 1) + ': ' + compStart + '–' + compEnd);
}

print('Composite Dataset Structure', labels.join('\n'));

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// ROBUSTNESS AND SENSITIVITY ASSESSMENT ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////// NUMBER OF OBSERVATIONS AND COVERAGE ///////////////////

//////////// NUMBER OF OBSERVATIONS ////////////////

// -----------------------
// 1) Per-pixel observation count for each multi-year interval
// -----------------------
function obsCountForInterval(startYear) {
  var start = ee.Date.fromYMD(startYear, 1, 1);
  var end = start.advance(STEP, 'year').advance(-1, 'day');
  var coll = landsat_collection_merged.filterDate(start, end);

  // Per-pixel count of valid scenes that went into the composite.
  // Any reflectance band works since all SR bands share the same mask after masksr().
  var count = coll.select('SR_B5')
                  .count()                      // counts images per pixel
                  .rename('N_obs')
                  .toInt16()
                  .set('start_year', startYear)
                  .set('end_year', end.get('year'))
                  .clip(roi);
  return count;
}

// Create an ImageCollection of observation-count images (one image per composite interval)
var obsCountsCollection = ee.ImageCollection(
  multiYearIntervals.map(obsCountForInterval)
);

// Observation count visualisation
// Adjusted max based on expected max observations per composite
// Approximate combined revisit interval by collection:
// dual-sensor combinations (offset orbits) ~8 days, single-sensor collections ~16 days
var REVISIT_DAYS = { 'A': 8, 'B': 8, 'C': 16, 'D': 16 };

// Theoretical maximum number of scenes = (summer days x STEP) / revisit interval
// Automatically adapts to STEP, POS_START/END, LANDSAT_COLLECTION
var computeVisMax = function() {
  var monthsInWindow = POS_END - POS_START + 1;
  var totalDays = monthsInWindow * 30.4 * STEP;
  var revisit = REVISIT_DAYS[LANDSAT_COLLECTION];
  return Math.ceil(totalDays / revisit);
};

var nVisMax = computeVisMax();

var nVis = {min: 0, max: nVisMax, palette: [
    '#d73027',  // 0–3   red    = very sparse
    '#fdae61',  // 4–9   orange = sparse
    '#fee08b',  // 10–15 yellow = moderate
    '#66bd63',  // 16–20 green  = good
    '#1a9850'   // 21+   dark green = dense
  ]
};

var meanN = obsCountsCollection.mean().rename('mean_N_obs');
Map.addLayer(meanN, {min: 0, max: nVisMax, palette: ['d73027','fee08b','1a9850']}, 'Mean N obs across all composites', false);

var list = obsCountsCollection.toList(obsCountsCollection.size());
for (var i = 0; i < N_COMPOSITES; i++) {
  var img = ee.Image(list.get(i));
  Map.addLayer(img, nVis, 'N obs ' + (START_YEAR + i * STEP), false);
}

// var stdN  = obsCountsCollection.reduce(ee.Reducer.stdDev()).rename('std_N');
// var cvN   = stdN.divide(meanN).multiply(100).rename('cv_N_obs');
// Map.addLayer(cvN, {min: 0, max: 100, palette: ['1a9850','fee08b','d73027']}, 'Temporal CV of N obs (%)', false);

// -----------------------
// 2) Total number of valid pixels (N_obs > 0) inside ROI for each composite
// -----------------------

var countsFC = ee.FeatureCollection(
  multiYearIntervals.map(function(y) {
    y = ee.Number(y);
    var img = obsCountForInterval(y);
    
    var pixelsWithObs = img.gt(0).reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: roi,
      scale: 30,            
      maxPixels: 1e13,
      tileScale: 4
    }).get('N_obs');

    return ee.Feature(null, {
      'start_year': y,
      'pixel_count': pixelsWithObs
    });
  })
);


// -----------------------
// 3) Build per-month observation-count for each multi-year interval
// -----------------------

// Build list of peak of season (POS) months in the aggregation window automatically

var POS_MONTHS = ee.List.sequence(POS_START, POS_END);
var POSMonths = POS_MONTHS.getInfo();

var monthBands = POSMonths.map(function(m) {
  return 'N_obs_' + m;
});

function monthlyObsForInterval(startYear) {
  startYear = ee.Number(startYear);
  var start = ee.Date.fromYMD(startYear, 1, 1);
  var end   = start.advance(STEP, 'year').advance(-1, 'day');
  var coll  = landsat_collection_merged.filterDate(start, end);
  
  var monthlyBands = POSMonths.map(function(m) {
    return coll.filter(ee.Filter.calendarRange(m, m, 'month'))
      .select('SR_B5')
      .count()
      .rename('N_obs_' + m);
  });

  var stacked = ee.ImageCollection(monthlyBands).toBands()
    .rename(monthBands)
    .set('start_year', startYear)
    .set('end_year',   end.get('year'))
    .clip(roi);

  return stacked;
}

// Create an ImageCollection of the monthly-observation images (one image per composite)
var monthlyObsCollection = ee.ImageCollection(multiYearIntervals.map(monthlyObsForInterval));
// print('Monthly obs per composite (ImageCollection):', monthlyObsCollection);

// -----------------------
// 4) Aggregate within ROI: count number of pixels (>=1 obs) per month per composite
// -----------------------

var monthlyCountsFC = ee.FeatureCollection(
  multiYearIntervals.map(function(y) {
    y = ee.Number(y);
    var img = monthlyObsForInterval(y);

    var binaryMonths = img.select(monthBands).gt(0);

    var counts = binaryMonths.reduceRegion({
      reducer:    ee.Reducer.sum(),
      geometry:   roi,
      scale:      30,
      maxPixels:  1e13,
      tileScale:  8
    });

    var props = { start_year: y };
    monthBands.forEach(function(bandName) {
      props['pixels_' + bandName] = counts.get(bandName);
    });

    return ee.Feature(null, props);
  })
);


// -----------------------
// Export the monthly summary (one row per composite)
// -----------------------

if (RUN_EXPORTS.monthlyCounts) {
  Export.table.toDrive({
    collection: monthlyCountsFC,   // FeatureCollection
    description: 'MonthlyPixelCounts_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION,
    folder: 'GEE_exports',                          
    fileNamePrefix: 'MonthlyPixelCounts_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION,
    fileFormat: 'CSV'
  });
}

// -----------------------
// Export the total valid-pixels-per-composite
// -----------------------

if (RUN_EXPORTS.totalCounts) {
  Export.table.toDrive({
    collection: countsFC,  // total pixels with any observation
    description: 'TotalValidPixelCount_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION,
    folder: 'GEE_exports',
    fileNamePrefix: 'TotalValidPixelCount_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION,
    fileFormat: 'CSV'
  });
}

//////////////////////////// COMPOSITE COMPLETENESS (Spatial Coverage) ///////////////

// -----------------------
// 1) Total pixel count in ROI: computed once, reused for every composite (INCLUDES WATER SURFACE)
// -----------------------

var totalPixelsInROI = ee.Image.constant(1).reduceRegion({
  reducer:   ee.Reducer.count(),
  geometry:  roi,
  scale:     30,
  maxPixels: 1e13,
  tileScale: 4
}).get('constant');

// -----------------------
// 2) Total pixel count in ROI: computed once, reused for every composite (EXCLUDES WATER SURFACE)
// -----------------------

// Land mask from JRC Global Surface Water — land = never observed as water.
// Matches the definition already used inside wtrmask() for consistency.
var landMask = GlobalWater.select('occurrence').gt(0).unmask(0).not().rename('land');

// Total land pixels in ROI — denominator for land-based completeness
var totalLandPixelsInROI = landMask.selfMask().reduceRegion({
  reducer:   ee.Reducer.count(),
  geometry:  roi,
  scale:     30,
  maxPixels: 1e13,
  tileScale: 4
}).get('land');


// -----------------------
// 3) Spatial coverage per composite
// -----------------------

function completenessForComposite(image) {
  
  var ndvi = image.select('NDVI');
  
  var toReduce = ee.Image([]);
  if (RUN_EXPORTS.totalCoverage) {
    toReduce = toReduce.addBands(ndvi.rename('valid_roi'));
  }
  if (RUN_EXPORTS.landCoverage) {
    toReduce = toReduce.addBands(ndvi.updateMask(landMask).rename('valid_land'));
  }
  
  var props = {
    'start_year': image.get('start_year'),
    'end_year':   image.get('end_year')
  };
  
  // Valid pixels anywhere in ROI
  var counts = toReduce.reduceRegion({
    reducer:   ee.Reducer.count(),
    geometry:  roi,
    scale:     30,
    maxPixels: 1e13,
    tileScale: 4
  })
  
  if (RUN_EXPORTS.totalCoverage) {
    var validCount = ee.Number(counts.get('valid_roi'));
    var total      = ee.Number(totalPixelsInROI);
    props.valid_pixels = validCount;
    props.total_pixels = total;
    props.completeness_pct = validCount.divide(total).multiply(100);
  }
  
  // Valid pixels restricted to land
  
  if (RUN_EXPORTS.landCoverage) {
    var validLandCount = ee.Number(counts.get('valid_land'));
    var totalLand      = ee.Number(totalLandPixelsInROI);
    props.valid_land_pixels = validLandCount;
    props.total_land_pixels = totalLand;
    props.land_completeness_pct = validLandCount.divide(totalLand).multiply(100);
  }

  return ee.Feature(null, props);
}

var completenessFC = ee.FeatureCollection(multiYearCollection.map(completenessForComposite));

// -----------------------
// 4) Export completeness stats
// -----------------------

if (RUN_EXPORTS.totalCoverage) {
  Export.table.toDrive({
    collection: completenessFC,                             
    description: 'TotalCompositeCompleteness_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION ,
    folder: 'GEE_exports',
    fileNamePrefix: 'TotalCompositeCompleteness_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION ,
    fileFormat: 'CSV',
    selectors:  ['start_year', 'end_year', 'total_pixels', 'valid_pixels', 'completeness_pct']
  });
}

if (RUN_EXPORTS.landCoverage) {
  Export.table.toDrive({
    collection: completenessFC,                             
    description: 'LandCompositeCompleteness_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION ,
    folder: 'GEE_exports',
    fileNamePrefix: 'LandCompositeCompleteness_' + 'POS' +  firstMonth + lastMonth + "_" + STEP + "yr_Colllection" + LANDSAT_COLLECTION ,
    fileFormat: 'CSV',
    selectors:  ['start_year', 'end_year', 'total_land_pixels', 'valid_land_pixels', 'land_completeness_pct']
  });
}

// -----------------------
// 5) Composites Visualisation
// -----------------------

// Loop through years and add layers dynamically (FOR VISUALIZATION)
var years = multiYearIntervals;
years.evaluate(function(yearsList){
  yearsList.forEach(function(year) {
    //var sampleroi = smallsample;
    var composite = multiYearCollection.filter(ee.Filter.eq('start_year', year)).first();
    // Map.addLayer(composite.select('NDVI'), {}, 'NDVI Composite ' + year, false);
    Map.addLayer(composite, nirVis, 'Composite ' + year, false);
    if (RUN_EXPORTS.compositeImages) {
      Export.image.toDrive({
        image: composite, //.select('band') //.select('NDVI')
        description: 'Collection' + LANDSAT_COLLECTION + '_' + STEP + 'yr_' + compstrat + '_Composite' + year + '_POS' + firstMonth + lastMonth, 
        scale: 30, 
        region: roi, 
        crs: 'EPSG:4326', 
        folder: 'GEE_Composites', 
        fileFormat: 'GeoTIFF',
        maxPixels: 1e13,
        formatOptions: {
        cloudOptimized: true
        }
      })
    }
  });
});


// Legend for N_obs
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px',
    backgroundColor: 'rgba(255,255,255,0.8)'
  }
});

var legendTitle = ui.Label({
  value: 'Observation count',
  style: {
    fontWeight: 'bold',
    fontSize: '14px',
    margin: '0 0 4px 0',
    padding: '0'
  }
});

legend.add(legendTitle);

// Helper to make a row
var makeRow = function(color, name) {
  var colorBox = ui.Label({
    style: {
      backgroundColor: '#' + color,
      padding: '8px',
      margin: '0 0 4px 0'
    }
  });

  var description = ui.Label({
    value: name,
    style: { margin: '0 0 4px 6px' }
  });

  return ui.Panel({
    widgets: [colorBox, description],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
};

var palette = ['d73027', 'fdae61', 'fee08b', '66bd63', '1a9850'];
var names = ['Low (0-3)','Sparse (4–6)', 'Moderate (7–9)', 'Good (10-12)', 'Dense (13+)']; 

for (var i = 0; i < palette.length; i++) {
  legend.add(makeRow(palette[i], names[i]));
}

Map.add(legend);

/// Map IFL feature collection for faster IFL ID identification

var inspector = ui.Label({
  value: 'Click a IFL to see its ID',
  style: {position: 'bottom-left', padding: '8px', fontWeight: 'bold'}
});
Map.add(inspector);

// 4. Set up the click callback event
Map.onClick(function(coords) {
  // Clear map of temporary selections if needed
  var point = ee.Geometry.Point(coords.lon, coords.lat);
  
  // Find the feature that intersects the click point
  var selected = NAM_IFL.filterBounds(point).first();
  
  // Use evaluate to pass the server-side property string to the client-side UI
  selected.get('IFL_ID').evaluate(function(result) {
    if (result) {
      inspector.setValue('IFL ID: ' + result);
    } else {
      inspector.setValue('No IFL selected');
    }
  });
});

Map.addLayer(NAM_IFL, {color: 'darkgreen'}, 'NAM IFL', false)

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// Compute NDVI Timeseries /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to compute mean NDVI for each image

var NDVIcollection = multiYearCollection.select('NDVI');

var meanNDVIlist = NDVIcollection.map(function(image) {
  var meanNDVI = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 30,
    maxPixels: 1e13
  }).get('NDVI');
  return ee.Feature(null, { 'year': image.get('start_year'), 'meanNDVI': meanNDVI });
});

// Convert to a list
// var meanNDVIvalues = meanNDVIlist.aggregate_array('meanNDVI');

if (RUN_EXPORTS.meanNDVI) {
  Export.table.toDrive({
    collection: meanNDVIlist,
    folder: 'GEE_Composites',
    description: 'MeanNDVI_' + compstrat + '_'  + sensors.getInfo() + '_' + STEP + 'yrComposites' + '_Summer' + POS_START + POS_END,
    fileFormat: 'CSV'
  });
}