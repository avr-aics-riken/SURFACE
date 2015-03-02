//
// Install
//   npm install ws
//   npm install socket.io
//

var WebSocketServer = require('ws').Server;
var http = require('http');
var app = http.createServer(handler)
var io = require('socket.io').listen(app, {'log level': 2})
var fs = require('fs');

var port = 8080;
var wsport = 8082;

function generateDataURI(mime, data)
{
  var datauri = 'data:' + mime + ';base64,' + data;
  return datauri;
}

var server = http.createServer();
var wss = new WebSocketServer({server: server, path: '/'});
wss.on('connection', function(ws) {
    console.log('connected');
    ws.on('message', function(data, flags) {
        if (flags.binary) {
          console.log("got binary");
          console.log("length: " + data.length);
          ws.send('ack');

          // Assume data = jpeg binary.
          var buffer = new Buffer(data, 'binary').toString('base64'); 
          var dataURI = generateDataURI("image/jpeg", buffer);
          //console.log(dataURI);
          io.emit("image", dataURI);

        } else {
          console.log('>>> ' + data);
          if (data == 'goodbye') { console.log('<<< galaxy'); ws.send('galaxy'); }
          if (data == 'hello') { console.log('<<< world'); ws.send('world'); }
        }
    });
    ws.on('close', function() {
      console.log('Connection closed!');
    });
    ws.on('error', function(e) {
    });
});
server.listen(wsport);
console.log('Listening websocket on  port ' + wsport + '...');

io.sockets.on('connection', function (socket) {

  console.log('server: sock io connect');

});

app.listen(port);

function handler (req, res) {
  //console.log(req)
  var filename = req.url
  var tmp = req.url.split('.');
  var type = tmp[tmp.length-1];

  if (req.url == '/') {
    filename = '/index.html'
  }
  fs.readFile(__dirname + filename,
  function (err, data) {
    if (err) {
      res.writeHead(500);
      return res.end('Error loading ' + filename);
    }
    switch (type) {
      case 'html':
        res.writeHead(200, {'Content-Type': 'text/html'});
        break;
      case 'js':
        //console.log('js:' + req.url)
        res.writeHead(200, {'Content-Type': 'text/javascript'});
        break;
      case 'css':
        res.writeHead(200, {'Content-Type': 'text/css'});
        break;
    }
    res.end(data);
  });
}

