# -*- coding: utf-8 -*-
# Store static resources

import pkg_resources
import base64


def getresource(path):
    return pkg_resources.resource_string(__name__, path).decode()


def imgtobase64(path):
    return base64.b64encode(pkg_resources.resource_string(__name__, path)).decode()


# HTML template
_html = {
    'page-top': getresource('static/html/page-top.html'),
    'bk-css': getresource('static/css/bk-css.css'),
    'page-bottom': getresource('static/html/page-bottom.html'),
    'network': getresource('static/html/network.html'),
    'species-top': getresource('static/html/species-top.html'),
    'species-each': getresource('static/html/species-each.html'),
    'species-bottom': getresource('static/html/species-bottom.html'),
    'reactions-top': getresource('static/html/reactions-top.html'),
    'reactions-each': getresource('static/html/reactions-each.html'),
    'reactions-bottom': getresource('static/html/reactions-bottom.html'),
    'narrowurl':  getresource('static/html/narrowurl.html'),
    'speciessvg-top': getresource('static/html/speciessvg-top.html'),
    'speciessvg-each': getresource('static/html/speciessvg-each.html'),
    'speciessvg-bottom': getresource('static/html/speciessvg-bottom.html'),
    'script-hidereac': getresource('static/js/script-hidereac.js'),
    'script-hidespec': getresource('static/js/script-hidespec.js'),
    'tr-top': getresource('static/html/tr-top.html'),
    'tr-specnone': getresource('static/html/tr-specnone.html'),
    'tr-reacnone': getresource('static/html/tr-reacnone.html'),
    'tr-bottom': getresource('static/html/tr-bottom.html'),
}


_static_js = {
    "jquery.min.js": getresource('static/js/jquery.min.js'),
    "bootstrap.bundle.min.js": getresource('static/js/bootstrap.bundle.min.js'),
    "jquery.easing.min.js": getresource('static/js/jquery.easing.min.js'),
    "scrollreveal.min.js": getresource('static/js/scrollreveal.min.js'),
    "jquery.magnific-popup.min.js": getresource('static/js/jquery.magnific-popup.min.js'),
    "creative.min.js": getresource('static/js/creative.min.js'),
}
_static_css = {
    "bootstrap.min.css": getresource('static/css/bootstrap.min.css'),
    "creative.min.css": getresource('static/css/creative.min.css'),
    "magnific-popup.min.css": getresource('static/css/magnific-popup.min.css'),
}
_static_img = {
    "fire.jpg": imgtobase64('static/img/fire.png'),
    'img-title': imgtobase64('static/img/img-title.png'),
}
