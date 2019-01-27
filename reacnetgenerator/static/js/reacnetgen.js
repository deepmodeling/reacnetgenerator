/* global jsnx */
/* global linkreac */
var canvas = $(document)[0].getElementById("canvas");

var G = new jsnx.Graph();
var timer = null;

jsnx.draw(G, {
    "element": canvas,
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
        "title" (d) { return d.label; },
        "xlink:href" (d) { return "#" + d.node + "_border"; },
        "onclick" (d) {
            return "clearTimeout(timer);timer = setTimeout(function(){addnode('" + d.node + "')}, 300);";
        },
        "ondblclick" (d) {
            return "clearTimeout(timer);G.removeNode('" + d.node + "');";
        },
        "width": 100,
        "height": 100,
        "x": -50,
        "y": -50,
    },
    "nodeStyle": {
        "border": "1px solid #ddd"
    },
    "edgeStyle": {
        "fill": "#999"
    },
    "edgeAttr": {
        "ondblclick" (d) {
            return "G.removeEdge('" + d.edge[0] + "','" + d.edge[1] + "');";
        },
    }
}, true);
$("#canvassec").addClass("mfp-hide");

$(".popup-modal").magnificPopup({
    "type": "inline",
    "preloader": false,
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
