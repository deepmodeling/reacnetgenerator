//CSS
import 'bootstrap/dist/css/bootstrap.min.css';
import 'magnific-popup/dist/magnific-popup.css';
import 'startbootstrap-creative/css/creative.min.css';
import './reacnetgen.css';

/* global linkreac */
global.$=global.jQuery=require('jquery');
require('bootstrap');
require('jquery.easing');
require('scrollreveal');
require("magnific-popup");
require('startbootstrap-creative/js/creative.min');
var jsnx = require("jsnetworkx");
var G = new jsnx.Graph();

$(function () {
    var canvas = $(document)[0].getElementById("canvas");
    jsnx.draw(G, {
        "element": canvas,
        d3: require("d3"),
        "layoutAttr": {
            "charge": -1000,
            "linkDistance": 300,
            "gravity": 0.05,
        },
        "panZoom": {
            "enabled": false
        },
        "nodeShape": "use",
        "nodeAttr": {
            "title"(d) { return d.label; },
            "xlink:href"(d) { return "#" + d.node + "_border"; },
            "ondblclick"(d) {
                return "clearTimeout(timer1);G.removeNode('" + d.node + "');";
            },
            "onmousedown": "isdrag = false;timer2 = setTimeout(function(){isdrag = true}, 300);",
            "onmouseup"(d) {
                return "if(!isdrag){clearTimeout(timer2);clearTimeout(timer1);timer1 = setTimeout(function(){addnode('" + d.node + "')}, 300);}else{isdrag=false}";
            },
            "width": 100,
            "height": 100,
            "x": -50,
            "y": -50,
        },
        "nodeStyle": {
            "border": "1px solid #ddd",
        },
        "edgeStyle": {
            "fill": "#999"
        },
        "edgeAttr": {
            "ondblclick"(d) {
                return "G.removeEdge('" + d.edge[0] + "','" + d.edge[1] + "');";
            },
        }
    }, true);
    $("#canvassec").addClass("mfp-hide");

    $(".popup-modal").magnificPopup({
        "type": "inline",
        "preloader": false,
    });
});

/**
* add nodes for the specified species
*/
function addnode(spec) {
    G.addNode(spec);
    if (spec in linkreac) {
        var rightspecs = linkreac[spec];
        for (var i = 0; i < rightspecs.length; i++) {
            G.addNode(rightspecs[i]);
            G.addEdge(rightspecs[i], spec);
        }
    }
}

//define global
window.$ = $
window.addnode = addnode;
window.G = G;
window.timer1 = null;
window.timer2 = null;
window.isdrag = false;
