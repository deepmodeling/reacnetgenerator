/*!
 * ReacNetGenerator (https://reacnetgenerator.njzjz.win/)
 * Copyright 2018-2019 East China Normal University
 */
//CSS
import 'bootstrap/dist/css/bootstrap.min.css';
import 'magnific-popup/dist/magnific-popup.css';
import 'startbootstrap-creative/css/creative.min.css';
import 'paginationjs/dist/pagination.css';
import 'bootstrap-select/dist/css/bootstrap-select.min.css';
import './reacnetgen.css';

/* global linkreac */
global.$ = global.jQuery = require('jquery');
require('bootstrap');
require('jquery.easing');
require('jsrender');
require('paginationjs');
require("magnific-popup");
require("bootstrap-select");
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
        "nodeShape": "image",
        "nodeAttr": {
            "title"(d) { return d.label; },
            "href"(d) { 
                var circle = '<circle cx="50" cy="50" r="45" stroke="#00f" stroke-width="2" fill="#fff" />';
                return "data:image/svg+xml;base64," + window.btoa('<svg class="spec" version="1.1" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">' + circle + speciessvg[d.node] + '</svg>');
            },
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
            "border": "1px solid #ddd"
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
    $("#canvasbutton").show();

    if (species.length > 1) {
        var timelist = [{ "value": 1, "text": "All" }]
        for (var i = 0; i < species.length; i++) {
            timelist.push({ "value": i + 2, "text": "Time " + String(i + 1) });
        }
        $("#timeselect").html($.templates("<option value={{:value}}>{{:text}}</option>").render(timelist));
        $("#timeselectli").removeClass("d-none");
        $("select#timeselect").change(function () {
            showresults($(this).val());
        });
    }
    if (reactionsabcd.length > 0) {
        $("#reactionsabcd").removeClass("d-none");
    }
    showresults(1);
});

/** 
 * show results
 */
function showresult(data, size, tmpl, result, pager) {
    $(pager).pagination({
        dataSource: data,
        pageSize: size,
        callback: function (data, pagination) {
            for (var i in data) {
                data[i]["svg"] = {};
                var p = ['s', 'l', 'r'];
                for (var ii in p) {
                    if (p[ii] in data[i]) {
                        if (typeof (data[i][p[ii]]) == "string") {
                            data[i]["svg"][data[i][p[ii]]] = speciessvg[data[i][p[ii]]];
                        } else {
                            for (var j in data[i][p[ii]]) {
                                data[i]["svg"][data[i][p[ii]][j]] = speciessvg[data[i][p[ii]][j]];
                            }
                        }

                    }
                }
            }
            $(result).html($.templates(tmpl).render(data));
            $(".popup-modal").magnificPopup({
                "type": "inline",
                "preloader": false,
            });
        }
    });
}

function showresults(time) {
    $("#networkresult").html(network[time - 1]);
    showresult(species[time - 1], speciesshownum, "#specTmpl", "#speciesresult", "#speciespager");
    showresult(reactions[time - 1], reactionsshownum, "#reacTmpl", "#reactionsresult", "#reactionspager");
    showresult(reactionsabcd, reactionsshownum, "#reacabcdTmpl", "#reactionsabcdresult", "#reactionsabcdpager");
    // select
    $("#speciesselect").html($.templates("<option value={{:s}}>{{:s}}</option>").render(species[time - 1]));
    $("#reactionsselect").html($.templates("<option value={{:s}}>{{:s}}</option>").render(species[time - 1]));
    $("#reactionsabcdselect").html($.templates("<option value={{:s}}>{{:s}}</option>").render(species[time - 1]));
    $("select#speciesselect").change(function () {
        if ($(this).val().length > 0) {
            var speciessearch = [];
            for (var i in species[time - 1]) {
                if ($(this).val().indexOf(species[time - 1][i]) >= 0) {
                    speciessearch.push(species[time - 1][i]);
                }
            }
        } else {
            speciessearch = species[time - 1];
        }
        showresult(speciessearch, speciesshownum, "#specTmpl", "#speciesresult", "#speciespager");
    });
    $("select#reactionsselect").change(function () {
        if ($(this).val().length > 0) {
            var reactionssearch = [];
            for (var i in reactions[time - 1]) {
                if ($(this).val().indexOf(reactions[time - 1][i]["l"][0]) >= 0 || $(this).val().indexOf(reactions[time - 1][i]["r"][0]) >= 0) {
                    reactionssearch.push(reactions[time - 1][i]);
                }
            }
        } else {
            reactionssearch = reactions[time - 1];
        }
        showresult(reactionssearch, reactionsshownum, "#reacTmpl", "#reactionsresult", "#reactionspager");
    });
    $("select#reactionsabcdselect").change(function () {
        if ($(this).val().length > 0) {
            var reactionsabcdsearch = [];
            for (var i in reactionsabcd) {
                var b = false;
                for (var j in reactionsabcd[i]['l']) {
                    if ($(this).val().indexOf(reactionsabcd[i]["l"][j]) >= 0) b = true;
                }
                for (var j in reactionsabcd[i]['r']) {
                    if ($(this).val().indexOf(reactionsabcd[i]["r"][j]) >= 0) b = true;
                }
                if (b) {
                    reactionsabcdsearch.push(reactionsabcd[i]);
                }
            }
        } else {
            reactionsabcdsearch = reactionsabcd;
        }
        showresult(reactionsabcdsearch, reactionsshownum, "#reacabcdTmpl", "#reactionsabcdresult", "#reactionsabcdpager");
    });
}

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

function savesvg() {
    var svgData = '<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">' + $(".jsnx")[0].outerHTML + "</svg>";
    var svgBlob = new Blob([svgData], { type: "image/svg+xml;charset=utf-8" });
    var svgUrl = URL.createObjectURL(svgBlob);
    var a = document.createElement('a');
    var filename = 'network.svg';
    a.href = svgUrl;
    a.download = filename;
    a.click();
    window.URL.revokeObjectURL(svgUrl);
}

function clearnode() {
    var nodes = G.nodes();
    for (var i = 0; i < nodes.length; i++) {
        G.removeNode(nodes[i]);
    }
}

//define global
window.$ = $
window.addnode = addnode;
window.savesvg = savesvg;
window.clearnode = clearnode;
window.G = G;
window.timer1 = null;
window.timer2 = null;
window.isdrag = false;
