// SPDX-License-Identifier: LGPL-3.0-or-later
/*!
 * ReacNetGenerator (https://reacnetgenerator.njzjz.win/)
 * Copyright 2018-2019 East China Normal University
 */

const {searchspecies, searchreaction} = require("./select.js");
const {getFormula} = require("./formula.js");

// CSS
/// #if process.env.REACNETGENERATOR_BUILDWEB
import './reacnetgen_web.scss'
/// #else
import './reacnetgen.scss'
/// #endif

/* global rngdata */
global.$ = global.jQuery = require('jquery');
global.regeneratorRuntime = require("regenerator-runtime");
window.bootstrap = require('bootstrap');
require('@popperjs/core');
global.anime = window.anime = require('animejs');
require('jsrender');
require('paginationjs');
require("magnific-popup");
require("bootstrap-select");
require("smiles-drawer");
/// #if !process.env.REACNETGENERATOR_BUILDWEB
require('startbootstrap-creative/dist/js/scripts');
/// #endif
var jsnx = require("@njzjz/jsnetworkx");
var jsnx = jsnx.default || jsnx;
var G = new jsnx.Graph();

$(function() {
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
    "element" : canvas,
    d3 : require("d3"),
    "layoutAttr" : {
      "charge" : -1000,
      "linkDistance" : 300,
      "gravity" : 0.05,
    },
    "panZoom" : {"enabled" : false},
    "nodeShape" : "image",
    "nodeAttr" : {
      "title"(d) { return d.label; },
      "xlink:href"(d) {
        var circle =
            '<circle cx="50" cy="50" r="45" stroke="#00f" stroke-width="2" fill="#fff" />';
        return "data:image/svg+xml;base64," +
               window.btoa(unescape(encodeURIComponent(
                   '<svg class="spec" version="1.1" viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">' +
                   circle + getSpecSvg(d.node) + '</svg>')));
      },
      "ondblclick"(
          d) { return `clearTimeout(timer1);G.removeNode('${d.node}');`; },
      "onmousedown" :
          "isdrag = false;timer2 = setTimeout(function(){isdrag = true}, 300);",
      "onmouseup"(d) {
        return `if(!isdrag){clearTimeout(timer2);clearTimeout(timer1);timer1 = setTimeout(function(){addnode('${
            d.node}')}, 300);}else{isdrag=false}`;
      },
      "width" : 100,
      "height" : 100,
      "x" : -50,
      "y" : -50,
    },
    "nodeStyle" : {"border" : "1px solid #ddd"},
    "edgeStyle" : {"fill" : "#999"},
    "edgeAttr" : {
      "ondblclick"(
          d) { return `G.removeEdge('${d.edge[0]}','${d.edge[1]}');`; },
    },
    "stickyDrag" : true,
  },
            true);
  $("#canvassec").addClass("mfp-hide");
  $("#canvasbutton").show();
}

function loadcitation() { $(".citation").html($("#citationTmpl").html()); }

function loadrngdata() {
  const text = $('#rngdata').html();
  if (handlejsondata(text)) {
    return;
  }
  // read from url
  const parsed = new URLSearchParams(location.search);
  const jdata = parsed.get("jdata");
  if (jdata) {
    $.get(decodeURIComponent(jdata), function(data) {
      if (!handlejsondata(data)) {
        addloadbutton();
      }
    }, 'text');
  } else {
    addloadbutton();
  }
}

// convert smiles to SVG
let smilesDrawer = new SmilesDrawer.SvgDrawer({
  height : 500,
  width : 500,
  bondThickness : 1,
  explicitHydrogens : true,
  compactDrawing : false,
});
let svgs = {};

const getSpecSvg = (smi) => svgs[smi];

const loadAllSpec =
    () => {
      $('.smiles[data-smiles]').each(function() {
        if (!$(obj).data("smiles-loaded")) {
          var obj = this;
          var spec = $(obj).data("smiles");
          storeSVG(spec, () => {
            var base64img =
                "data:image/svg+xml;base64," +
                window.btoa(unescape(encodeURIComponent(getSpecSvg(spec))));
            $(obj).html("<img src='" + base64img + "'/>");
            $(obj).data("smiles-loaded", true);
          });
        }
      });
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
    // load time select
    var timelist = [ {"value" : 1, "text" : "All"} ];
    timelist = timelist.concat([...rngdata['species'].keys() ].slice(1).map(
        ii => ({"value" : ii + 1, "text" : `Time ${ii}`})));
    $("#timeselect").html($.templates("#optionTimeTmpl").render(timelist));
    $("#timeselectli").removeClass("d-none");
    $("select#timeselect")
        .on("change", function() { showresults($(this).val()); });
  }
  if (rngdata['reactionsabcd'].length) {
    $("#reactionsabcd").removeClass("d-none");
  }
  showresults(1);
}

function loadsection() {
  // show sections
  const sections =
      [ 'network', 'species', 'reactions', 'reactionsabcd' ].filter(
          ii => rngdata[ii]);
  sections.forEach(ii => $(`#${ii}`).show());
  $("#navs").append($.templates("#navTmpl").render(sections));
  $("#buttons").html($.templates("#buttonTmpl").render(sections));

  // set anime again after new button appears
  $('a.js-scroll-trigger[href*="#"]:not([href="#"])').on('click', function() {
    if (location.pathname.replace(/^\//, "") ==
            this.pathname.replace(/^\//, "") &&
        location.hostname == this.hostname) {
      var target = $(this.hash);
      target = target.length ? target : $(`[name=${this.hash.slice(1)}]`);
      if (target.length) {
        anime({
          targets : 'html, body',
          scrollTop : target.offset().top - 72,
          duration : 1000,
          easing : 'easeInOutExpo'
        });
        return false;
      }
    }
  });
}

/**
 * show results
 */
function showresult(data, size, tmpl, result, pager) {
  $(pager).pagination({
    dataSource : data,
    pageSize : size,
    callback : function(data, pagination) {
      $(result).html($.templates(tmpl).render(data));
      loadAllSpec();
      $(".popup-modal").magnificPopup({
        "type" : "inline",
        "preloader" : false,
      });
    }
  });
}

function showresults(time) {
  const specdata = rngdata['species'][time - 1];
  const reactionsdata = rngdata['reactions'][time - 1];
  const reactionsabcddata = rngdata['reactionsabcd'];
  const speciesshownum = rngdata['speciesshownum'];
  const reactionsshownum = rngdata['reactionsshownum'];
  $("#networkresult").html(rngdata['network'][time - 1]);
  showresult(specdata, speciesshownum, "#specTmpl", "#speciesresult",
             "#speciespager");
  showresult(reactionsdata, reactionsshownum, "#reacTmpl", "#reactionsresult",
             "#reactionspager");
  showresult(reactionsabcddata, reactionsshownum, "#reacabcdTmpl",
             "#reactionsabcdresult", "#reactionsabcdpager");
  // select
  // render formula to show formula in the select
  specdata.forEach(dd => { dd['formula'] = getFormula(dd['s']); });
  // limit size of options to 65536
  // workaround to fix
  // https://github.com/snapappointments/bootstrap-select/issues/2793
  const specdata_minify = specdata.slice(0, 65536);
  $("#speciesselect").html($.templates("#optionTmpl").render(specdata_minify));
  $("#reactionsselect")
      .html($.templates("#optionTmpl").render(specdata_minify));
  $("#reactionsabcdselect")
      .html($.templates("#optionTmpl").render(specdata_minify));
  $("select#speciesselect").on("change", function() {
    const speciessearch = searchspecies($(this).val(), specdata);
    showresult(speciessearch, speciesshownum, "#specTmpl", "#speciesresult",
               "#speciespager");
  });
  $("select#reactionsselect").on("change", function() {
    const reactionssearch = searchreaction($(this).val(), reactionsdata);
    showresult(reactionssearch, reactionsshownum, "#reacTmpl",
               "#reactionsresult", "#reactionspager");
  });
  $("select#reactionsabcdselect").on("change", function() {
    const reactionsabcdsearch =
        searchreaction($(this).val(), reactionsabcddata);
    showresult(reactionsabcdsearch, reactionsshownum, "#reacabcdTmpl",
               "#reactionsabcdresult", "#reactionsabcdpager");
  });
  // refresh select picker
  $('.selectpicker').selectpicker("refresh");
}

/**
 * add nodes for the specified species
 */
function addnode(spec) {
  addSingleNode(spec);
  if (spec in rngdata['linkreac']) {
    rngdata['linkreac'][spec].forEach(rightspec => {
      addSingleNode(rightspec);
      G.addEdge(rightspec, spec);
    });
  }
}

function storeSVG(spec, callback) {
  // load smiles svg
  if (spec in svgs) {
    callback();
  } else {
    SmilesDrawer.parse(
        spec,
        function(tree) {
          smilesDrawer.draw(tree, 'tmpsvg', 'light', false);
          svgs[spec] = $('#tmpsvgdiv').html();
          callback();
        },
        function(err) {
          console.log(err);
          svgs[spec] =
              `<svg viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg"><text x="0" y="50">${
                  spec}</text></svg>`;
          callback();
        });
  }
}

function addSingleNode(spec) { storeSVG(spec, () => { G.addNode(spec); }); }

function savesvg() {
  var svgData =
      '<svg version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">' +
      $(".jsnx")[0].outerHTML + "</svg>";
  var svgBlob = new Blob([ svgData ], {type : "image/svg+xml;charset=utf-8"});
  var svgUrl = URL.createObjectURL(svgBlob);
  var a = document.createElement('a');
  var filename = 'network.svg';
  a.href = svgUrl;
  a.download = filename;
  a.click();
  window.URL.revokeObjectURL(svgUrl);
}

function clearnode() { G.nodes().foreach(node => G.removeNode(node)); }

function addloadbutton() {
  $("#buttons").html($("#loadTmpl").html());
  $('#loadbutton').on("change", function(e) {
    const f = e.target.files[0];
    const reader = new FileReader();
    reader.onload = (function(
        theFile) { return function(e) { handlejsondata(e.target.result); }; })(
        f);
    reader.readAsText(f);
  });
}

// placeholder for SimpleLightbox
function SimpleLightbox(config) {};

// define global
window.$ = $;
window.addnode = addnode;
window.savesvg = savesvg;
window.clearnode = clearnode;
window.G = G;
window.timer1 = null;
window.timer2 = null;
window.isdrag = false;
window.SimpleLightbox = SimpleLightbox;
