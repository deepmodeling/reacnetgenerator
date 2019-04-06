//CSS
import 'bootstrap/dist/css/bootstrap.min.css';
import 'magnific-popup/dist/magnific-popup.css';
import 'startbootstrap-creative/css/creative.min.css';
import './reacnetgen.css';

/* global linkreac */
global.$ = global.jQuery = require('jquery');
require('bootstrap');
require('jquery.easing');
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
        },
        "stickyDrag": true,
    }, true);
    $("#canvassec").addClass("mfp-hide");

    $(".popup-modal").magnificPopup({
        "type": "inline",
        "preloader": false,
    });
    $("select#timeselect").change(function(){
        $(".time").hide();
        $(".time-"+$(this).val()).show();
    });
    $(".time").hide();
    $(".time-1").show();
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

function savesvg(){
    var svgData = '<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">'+$(".jsnx")[0].outerHTML+$("#svgdefs")[0].outerHTML+$("#svgspecs")[0].outerHTML+"</svg>";
    var svgBlob = new Blob([svgData], {type:"image/svg+xml;charset=utf-8"});
    var svgUrl = URL.createObjectURL(svgBlob);
    var a = document.createElement('a');
    var filename = 'network.svg';
    a.href = svgUrl;
    a.download = filename;
    a.click();
    window.URL.revokeObjectURL(svgUrl);
}

//define global
window.$ = $
window.addnode = addnode;
window.savesvg = savesvg;
window.G = G;
window.timer1 = null;
window.timer2 = null;
window.isdrag = false;
