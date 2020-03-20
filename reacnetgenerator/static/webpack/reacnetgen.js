/*!
 * ReacNetGenerator (https://reacnetgenerator.njzjz.win/)
 * Copyright 2018-2019 East China Normal University
 */
//CSS
import './reacnetgen.scss'

/* global rngdata */
global.$ = global.jQuery = require('jquery');
global.regeneratorRuntime = require("regenerator-runtime");
require('bootstrap');
require('jquery.easing');
require('jsrender');
require('paginationjs');
require("magnific-popup");
require("bootstrap-select");
require('startbootstrap-creative/dist/js/scripts');
var jsnx = require("@njzjz/jsnetworkx");
var G = new jsnx.Graph();

$(function () {
    loadcitation();
    drawcanvas();
    loadrngdata();
});

function handlerngdata(rngdata) {
    global.rngdata = rngdata;
    loadsection();
    loaddata();
}

function drawcanvas() {
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
            "xlink:href"(d) {
                var circle = '<circle cx="50" cy="50" r="45" stroke="#00f" stroke-width="2" fill="#fff" />';
                return "data:image/svg+xml;base64," + window.btoa('<svg class="spec" version="1.1" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">' + circle + rngdata['speciessvg'][d.node] + '</svg>');
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
}

function loadcitation() {
    $(".citation").html($("#citationTmpl").html());
}

function loadrngdata() {
    var text = $('#rngdata').html();
    if (handlejsondata(text)) {
        return;
    }
    // read from url
	const queryString = require('query-string');
	const parsed = queryString.parse(location.search);
	var jdata = parsed['jdata']
    if (jdata) {
        $.get(decodeURIComponent(jdata), function (data) {
            if (!handlejsondata(data)) {
                addloadbutton();
            }
        }, 'text');
    } else {
        addloadbutton();
    }
}

function handlejsondata(text) {
    try {
        handlerngdata(JSON.parse(text));
        return true;
    } catch (err) {
        console.log(err);
        return false;
    }
}

function loaddata() {
    if (rngdata['species'].length > 1) {
        var timelist = [{ "value": 1, "text": "All" }]
        for (var i = 1; i < rngdata['species'].length; i++) {
            timelist.push({ "value": i + 1, "text": "Time " + String(i) });
        }
        $("#timeselect").html($.templates("#optionTimeTmpl").render(timelist));
        $("#timeselectli").removeClass("d-none");
        $("select#timeselect").change(function () {
            showresults($(this).val());
        });
    }
    if (rngdata['reactionsabcd'].length > 0) {
        $("#reactionsabcd").removeClass("d-none");
    }
    showresults(1);
}

function loadsection() {
    var sections = [];
    if (rngdata['network']) {
        sections.push('network');
        $('#network').show();
    }
    if (rngdata['species']) {
        sections.push('species');
        $('#species').show();
    }
    if (rngdata['reactions']) {
        sections.push('reactions');
        $('#reactions').show();
    }
    if (rngdata['reactionsabcd']) {
        $('#reactionsabcd').show();
    }
    $("#navs").append($.templates("#navTmpl").render(sections));
    $("#buttons").html($.templates("#buttonTmpl").render(sections));
}

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
                            data[i]["svg"][data[i][p[ii]]] = rngdata['speciessvg'][data[i][p[ii]]];
                        } else {
                            for (var j in data[i][p[ii]]) {
                                data[i]["svg"][data[i][p[ii]][j]] = rngdata['speciessvg'][data[i][p[ii]][j]];
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
    $("#networkresult").html(rngdata['network'][time - 1]);
    showresult(rngdata['species'][time - 1], rngdata['speciesshownum'], "#specTmpl", "#speciesresult", "#speciespager");
    showresult(rngdata['reactions'][time - 1], rngdata['reactionsshownum'], "#reacTmpl", "#reactionsresult", "#reactionspager");
    showresult(rngdata['reactionsabcd'], rngdata['reactionsshownum'], "#reacabcdTmpl", "#reactionsabcdresult", "#reactionsabcdpager");
    // select
    $("#speciesselect").html($.templates("#optionTmpl").render(rngdata['species'][time - 1]));
    $("#reactionsselect").html($.templates("#optionTmpl").render(rngdata['species'][time - 1]));
    $("#reactionsabcdselect").html($.templates("#optionTmpl").render(rngdata['species'][time - 1]));
    $("select#speciesselect").change(function () {
        if ($(this).val().length > 0) {
            var speciessearch = [];
            for (var i in rngdata['species'][time - 1]) {
                if ($(this).val().indexOf(rngdata['species'][time - 1][i]['s']) >= 0) {
                    speciessearch.push(rngdata['species'][time - 1][i]);
                }
            }
        } else {
            speciessearch = rngdata['species'][time - 1];
        }
        showresult(speciessearch, rngdata['speciesshownum'], "#specTmpl", "#speciesresult", "#speciespager");
    });
    $("select#reactionsselect").change(function () {
        if ($(this).val().length > 0) {
            var reactionssearch = [];
            for (var i in rngdata['reactions'][time - 1]) {
                if ($(this).val().indexOf(rngdata['reactions'][time - 1][i]["l"][0]) >= 0 || $(this).val().indexOf(rngdata['reactions'][time - 1][i]["r"][0]) >= 0) {
                    reactionssearch.push(rngdata['reactions'][time - 1][i]);
                }
            }
        } else {
            reactionssearch = rngdata['reactions'][time - 1];
        }
        showresult(reactionssearch, rngdata['reactionsshownum'], "#reacTmpl", "#reactionsresult", "#reactionspager");
    });
    $("select#reactionsabcdselect").change(function () {
        if ($(this).val().length > 0) {
            var reactionsabcdsearch = [];
            for (var i in rngdata['reactionsabcd']) {
                var b = false;
                for (var j in rngdata['reactionsabcd'][i]['l']) {
                    if ($(this).val().indexOf(rngdata['reactionsabcd'][i]["l"][j]) >= 0) b = true;
                }
                for (var j in rngdata['reactionsabcd'][i]['r']) {
                    if ($(this).val().indexOf(rngdata['reactionsabcd'][i]["r"][j]) >= 0) b = true;
                }
                if (b) {
                    reactionsabcdsearch.push(rngdata['reactionsabcd'][i]);
                }
            }
        } else {
            reactionsabcdsearch = rngdata['reactionsabcd'];
        }
        showresult(reactionsabcdsearch, rngdata['reactionsshownum'], "#reacabcdTmpl", "#reactionsabcdresult", "#reactionsabcdpager");
    });
    // refresh select picker
    $('.selectpicker').selectpicker("refresh");
}

/**
* add nodes for the specified species
*/
function addnode(spec) {
    G.addNode(spec);
    if (spec in rngdata['linkreac']) {
        var rightspecs = rngdata['linkreac'][spec];
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

function addloadbutton() {
    $("#buttons").html($("#loadTmpl").html());
    $('#loadbutton').change(function (e) {
        var f = e.target.files[0];
        var reader = new FileReader();
        reader.onload = (function (theFile) {
            return function (e) {
                handlejsondata(e.target.result);
            };
        })(f);
        reader.readAsText(f);
    });
}

//define global
window.$ = $;
window.addnode = addnode;
window.savesvg = savesvg;
window.clearnode = clearnode;
window.G = G;
window.timer1 = null;
window.timer2 = null;
window.isdrag = false;
