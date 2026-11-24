#!/usr/bin/env python3
"""A logging forwarding proxy for the DASH REST server: dash_client talks to http://127.0.0.1:PORT/dash/REST1.0/, the proxy forwards each request to the real server (https://sysimm.org/dash/REST1.0/) and appends one JSON
line per exchange to the log (method, path+query, request body and response body as base64, status, content-type), so that the next session can replay the server offline.  Standard library only; the
User-Agent towards the server is the neutral "Seqeron-review/1.0".  Usage: python3 dashproxy.py PORT LOGFILE [UPSTREAM]"""
import base64, json, sys, threading, urllib.request, urllib.error
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

PORT = int(sys.argv[1])
LOG = sys.argv[2]
UPSTREAM = (sys.argv[3] if len(sys.argv) > 3 else 'https://sysimm.org').rstrip('/')
LOCK = threading.Lock()
SEQ = [0]


class H(BaseHTTPRequestHandler):
    protocol_version = 'HTTP/1.0'

    def log_message(self, *a):
        pass

    def _go(self):
        n = int(self.headers.get('Content-Length') or 0)
        body = self.rfile.read(n) if n else b''
        url = UPSTREAM + self.path
        req = urllib.request.Request(url, data=body if self.command == 'POST' else None, method=self.command)
        req.add_header('User-Agent', 'Seqeron-review/1.0')
        if self.headers.get('Content-Type'):
            req.add_header('Content-Type', self.headers['Content-Type'])
        status, rbody, ctype = 0, b'', ''
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                status, rbody, ctype = r.status, r.read(), r.headers.get('Content-Type', '')
        except urllib.error.HTTPError as e:
            status, rbody, ctype = e.code, e.read(), e.headers.get('Content-Type', '')
        except Exception as e:  # noqa
            status, rbody, ctype = 599, ('proxy error: %s' % e).encode(), 'text/plain'
        with LOCK:
            SEQ[0] += 1
            with open(LOG, 'a') as f:
                f.write(json.dumps({'n': SEQ[0], 'method': self.command, 'path': self.path, 'request_b64': base64.b64encode(body).decode(), 'status': status,
                                    'content_type': ctype, 'response_b64': base64.b64encode(rbody).decode()}) + '\n')
        try:  # the client may have closed the connection already (its own timeout); the exchange is logged above
            self.send_response(status)
            self.send_header('Content-Type', ctype or 'application/octet-stream')
            self.send_header('Content-Length', str(len(rbody)))
            self.end_headers()
            self.wfile.write(rbody)
        except (BrokenPipeError, ConnectionResetError):
            pass

    do_GET = _go
    do_POST = _go


if __name__ == '__main__':
    ThreadingHTTPServer(('127.0.0.1', PORT), H).serve_forever()
