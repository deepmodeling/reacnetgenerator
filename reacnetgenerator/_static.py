# -*- coding: utf-8 -*-
"""Store static resources."""

import base64

import pkg_resources

def _getresource(path):
    return pkg_resources.resource_string(__name__, path).decode()


def _imgtobase64(path):
    return base64.b64encode(
        pkg_resources.resource_string(__name__, path)).decode()

# HTML template
_html = {
    'template': _getresource('static/template.html'),
}

_static_js = {
    "jsnetworkx.js": _getresource('static/webpack/bundle.js'),
}

_static_img = {
    'img-title': _imgtobase64('static/img-title.png'),
}
