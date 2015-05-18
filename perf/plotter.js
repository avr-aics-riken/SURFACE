//
// Assume benchmark was ran with --color_print=false
//
require('shelljs/global');

var swig  = require('swig');
var fs = require('fs');

if (process.argv.length < 3) {
  console.log("needs input.dat");
  process.exit(1);
}

// -- config ------------------------------------------
var phantomjs = '../tools/macosx/phantomjs'
// ----------------------------------------------------

// http://stackoverflow.com/questions/10645994/node-js-how-to-format-a-date-string-in-utc
function dateFormat (date, fstr, utc) {
  utc = utc ? 'getUTC' : 'get';
  return fstr.replace (/%[YmdHMS]/g, function (m) {
    switch (m) {
    case '%Y': return date[utc + 'FullYear'] (); // no leading zeros required
    case '%m': m = 1 + date[utc + 'Month'] (); break;
    case '%d': m = date[utc + 'Date'] (); break;
    case '%H': m = date[utc + 'Hours'] (); break;
    case '%M': m = date[utc + 'Minutes'] (); break;
    case '%S': m = date[utc + 'Seconds'] (); break;
    default: return m.slice (1); // unknown code, remove %
    }
    // add leading zero if required
    return ('0' + m).slice (-2);
  });
}

// 1024k -> 1M -> 1.0
// 2M    ->       2.0
function convertValue (elem) {

  if (elem.match(/k$/)) {
    var val = parseFloat(elem.slice(0, -1));
    var megaVal = (val / 1024.0)
    return megaVal;
  } else if (elem.match(/M$/)) {
    var val = parseFloat(elem.slice(0, -1));
    return val;
  } else {
    // ??? @todo
    return elem;
  }

}

var data = fs.readFileSync(process.argv[2]).toString();

var dict = {}

data.split('\n').forEach(function (line) {
  var re = line.match(/DEBUG:\s*(\S+)\/(\S+)\s*(\d+)\s*(\d+)/);
  if (re && re.length == 5) {
    var suite = re[1];
    var variant = re[2];
    var cpu_time = re[3] / (1000 * 1000 * 1000); // [ns] -> [sec]

    if (dict[suite] == undefined) {
      dict[suite] = []
    }

    dict[suite].push({elems: convertValue(variant), secs: cpu_time});

  }
});


for (var key in dict) {

  var suite = {}
  suite['suite'] = key;
  suite['data']  = dict[key]; // JSON

  var result = swig.renderFile(__dirname + '/template.html', suite);

  // generate filename 
  // Assume 'suite' string is safe to use as filename. 
  var date = new Date();
  var dateString = dateFormat(new Date(), "%Y-%m-%d", false);

  var htmlFilename = key + '-' + dateString + '.html';
  var pdfFilename = key + '-' + dateString + '.pdf';

  fs.writeFileSync(htmlFilename, result);
  console.log("Wrote  : " + htmlFilename);

  // Invoke phantom.js
  cmd = phantomjs + ' rasterize.js ' + htmlFilename + ' ' + pdfFilename + ' A4'
  if (exec(cmd).code != 0) {
    echo('Error: html -> pdf conversion using phantom.js failed: ' + cmd);
    exit(1);
  }
  console.log("Convert: " + pdfFilename);
};
