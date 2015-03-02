var fs = require('fs');
var exec = require('child_process').exec;

var lsglsandbox = './lsglsandbox'
var workdir = './'
var shadername = process.argv[3] || "input.frag"

fs.watchFile(shadername, function(curr, prev) {
  console.log('watch')

  var child = exec(lsglsandbox + " " + shadername, {cwd: workdir}, function (error, stdout, stderr) {
      console.log('stdout: ' + stdout);
      console.log('stderr: ' + stderr);
      if (error !== null) {
        console.log('exec error: ' + error);
      }
  });
  
}
