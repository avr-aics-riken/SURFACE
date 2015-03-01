//
// npm install ws
//

var WebSocketServer = require('ws').Server;
var http = require('http');

var port = 8080;

var server = http.createServer();
var wss = new WebSocketServer({server: server, path: '/'});
wss.on('connection', function(ws) {
    console.log('connected');
    ws.on('message', function(data, flags) {
        if (flags.binary) {
          console.log("got binary");
          console.log("length: " + data.length);
          ws.send('ack');
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
server.listen(port);
console.log('Listening websocket on  port ' + port + '...');

