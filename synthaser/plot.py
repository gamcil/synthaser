import http.server
import socketserver
import webbrowser
import json
import shutil
import pathlib
import logging

from functools import partial

from synthaser import grouping


LOG = logging.getLogger(__name__)


def get_data(container):
    container.sort(key=lambda s: (s.classification, -s.sequence_length))
    hierarchy = grouping.get_classification_paths(container)
    return {
        "synthases": {s.header: s.to_dict() for s in container},
        "order": [s.header for s in container],
        "types": dict(grouping.group_synthases(container)),
        "groups": grouping.get_annotation_groups(hierarchy)
    }


class CustomHandler(http.server.BaseHTTPRequestHandler):
    """Handler for serving cblaster plots."""

    def __init__(self, data, *args, **kwargs):
        self._data = data
        self._dir = pathlib.Path(__file__).resolve().parent.parent / "synthaser" / "plot"
        super().__init__(*args, **kwargs)

    def copy_file(self, source):
        shutil.copyfileobj(source, self.wfile)

    def send_headers(self, mime):
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.end_headers()

    def log_message(self, format, *args):
        """Suppresses logging messages on every request."""
        return

    def do_GET(self):
        """Serves each component of the cblaster plot."""
        if self.path == "/data.json":
            self.send_headers("text/json")
            self.wfile.write(json.dumps(self._data).encode())
            return
        path, mime = None, None
        if self.path == "/":
            path, mime = self._dir / "index.html", "text/html"
        elif self.path == "/index.css":
            path, mime = self._dir / "index.css", "text/css"
        elif self.path == "/d3.min.js":
            path, mime = self._dir / "d3.min.js", "text/javascript"
        elif self.path == "/synthaser.js":
            path, mime = self._dir / "synthaser.js", "text/javascript"
        elif self.path == "/synthaser.min.js":
            path, mime = self._dir / "synthaser.min.js", "text/javascript"
        if not path:
            return
        with path.open("rb") as fp:
            self.send_headers(mime)
            self.copy_file(fp)


def serve_html(data):
    """Serve a synthaser plot using the socketserver module."""
    handler = partial(CustomHandler, data)

    # Instantiate a new server, bind to any open port
    with socketserver.TCPServer(("localhost", 0), handler) as httpd:

        # Automatically open web browser to bound address
        address, port = httpd.server_address
        url = f"http://{address}:{port}/"
        webbrowser.open(url)

        # Start serving the plot; shutdown on a keyboard interrupt
        try:
            LOG.info(f"Serving synthaser plot at {url} (Ctrl+C to stop).")
            httpd.serve_forever()
        except KeyboardInterrupt:
            httpd.shutdown()


def save_html(data, output):
    """Generates a static HTML file with all visualisation code."""

    directory = pathlib.Path(__file__).resolve().parent.parent / "synthaser" / "plot"

    with (directory / "index.html").open() as fp:
        html = fp.read()

    css_string = '<link href="index.css" rel="stylesheet"></link>'
    d3_string = '<script src="d3.min.js"></script>'
    sy_string = '<script src="synthaser.js"></script>'
    sy_min = '<script src="synthaser.min.js"></script>'

    with (directory / "index.css").open() as fp:
        css = fp.read()
        html = html.replace(css_string, f"<style>{css}</style>")

    with (directory / "d3.min.js").open() as fp:
        d3 = fp.read()
        html = html.replace(d3_string, f"<script>{d3}</script>")

    with (directory / "synthaser.js").open() as fp:
        sy = f"const data={json.dumps(data)};" + fp.read()
        html = html.replace(sy_string, f"<script>{sy}</script>")

    with (directory / "synthaser.min.js").open() as fp:
        text = fp.read()
        html = html.replace(sy_min, f"<script>{text}</script>")

    with open(output, "w") as fp:
        fp.write(html)


def plot_synthases(container, output=None):
    """Generates synthaser plot from a collection of Synthase objects."""
    data = get_data(container)
    if output:
        save_html(data, output)
        webbrowser.open(output)
    else:
        serve_html(data)
