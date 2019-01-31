# -*- coding: utf-8 -*-
"""Store static resources."""

import base64

import pkg_resources
from cssmin import cssmin
from jsmin import jsmin


def _getresource(path):
    return pkg_resources.resource_string(__name__, path).decode()


def _imgtobase64(path):
    return base64.b64encode(
        pkg_resources.resource_string(__name__, path)).decode()


# HTML template
_html = {
    'bk-css': cssmin(_getresource('static/css/bk-css.css')),
    "reacnetgen.js": jsmin(_getresource('static/js/reacnetgen.js')),
    'template': _getresource('static/html/template.html'),
}

_static_js = {
    "jquery.min.js": _getresource('static/js/vendor/jquery.min.js'),
    "bootstrap.bundle.min.js": _getresource('static/js/vendor/bootstrap.bundle.min.js'),
    "jquery.easing.min.js": _getresource('static/js/vendor/jquery.easing.min.js'),
    "scrollreveal.min.js": _getresource('static/js/vendor/scrollreveal.min.js'),
    "jquery.magnific-popup.min.js": _getresource('static/js/vendor/jquery.magnific-popup.min.js'),
    "creative.min.js": _getresource('static/js/vendor/creative.min.js'),
    "d3.min.js": _getresource('static/js/vendor/d3.min.js'),
    "jsnetworkx.js": _getresource('static/js/vendor/jsnetworkx.min.js'),
}
_static_css = {
    "bootstrap.min.css": _getresource('static/css/vendor/bootstrap.min.css'),
    "creative.min.css": _getresource('static/css/vendor/creative.min.css'),
    "magnific-popup.min.css": _getresource('static/css/vendor/magnific-popup.min.css'),
}
_static_img = {
    "fire.jpg": _imgtobase64('static/img/fire.png'),
    'img-title': _imgtobase64('static/img/img-title.png'),
}
