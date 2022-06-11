(window.webpackJsonp=window.webpackJsonp||[]).push([[10],{438:function(a,t,e){"use strict";e.r(t);var n=e(65),s=Object(n.a)({},(function(){var a=this,t=a.$createElement,e=a._self._c||t;return e("ContentSlotsDistributor",{attrs:{"slot-key":a.$parent.slotKey}},[e("h1",{attrs:{id:"guide"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#guide"}},[a._v("#")]),a._v(" Guide")]),a._v(" "),e("p",[a._v("ReacNetGenerator is an automatic reaction network generator for reactive molecular dynamics simulation. This page provides a easy way to start it.")]),a._v(" "),e("h2",{attrs:{id:"installation"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#installation"}},[a._v("#")]),a._v(" Installation")]),a._v(" "),e("p",[a._v("You can "),e("a",{attrs:{href:"https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html",target:"_blank",rel:"noopener noreferrer"}},[a._v("install Anaconda or Miniconda"),e("OutboundLink")],1),a._v(" to obtain conda, and install ReacNetGenerator easily with conda:")]),a._v(" "),e("div",{staticClass:"language-bash extra-class"},[e("pre",{pre:!0,attrs:{class:"language-bash"}},[e("code",[a._v("conda "),e("span",{pre:!0,attrs:{class:"token function"}},[a._v("install")]),a._v(" reacnetgenerator -c conda-forge\n")])])]),e("p",[a._v("See "),e("RouterLink",{attrs:{to:"/guide/build.html"}},[a._v("the build guide")]),a._v(" if you want to build ReacNetGenerator by yourself.")],1),a._v(" "),e("h2",{attrs:{id:"usage"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#usage"}},[a._v("#")]),a._v(" Usage")]),a._v(" "),e("h3",{attrs:{id:"command-line"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#command-line"}},[a._v("#")]),a._v(" Command line")]),a._v(" "),e("p",[a._v("ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:")]),a._v(" "),e("div",{staticClass:"language-bash extra-class"},[e("pre",{pre:!0,attrs:{class:"language-bash"}},[e("code",[a._v("reacnetgenerator --type lammpsdumpfile -i dump.reaxc -a C H O --nohmm\n")])])]),e("p",[a._v("where C, H, and O are atomic names in the input file. "),e("a",{attrs:{href:"/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json",target:"_blank"}},[a._v("Analysis report")]),a._v(" will be generated automatically.")]),a._v(" "),e("p",[a._v("Also, ReacNetGenerator can process files containing bond information, e.g. LAMMPS bond file:")]),a._v(" "),e("div",{staticClass:"language-bash extra-class"},[e("pre",{pre:!0,attrs:{class:"language-bash"}},[e("code",[a._v("reacnetgenerator --type lammpsbondfile -i bonds.reaxc -a C H O --nohmm\n")])])]),e("p",[a._v("You can running the following script for help:")]),a._v(" "),e("div",{staticClass:"language-bash extra-class"},[e("pre",{pre:!0,attrs:{class:"language-bash"}},[e("code",[a._v("reacnetgenerator -h\n")])])]),e("h3",{attrs:{id:"enable-hmm-filter"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#enable-hmm-filter"}},[a._v("#")]),a._v(" Enable HMM filter")]),a._v(" "),e("p",[a._v("One can remove "),e("code",[a._v("--nohmm")]),a._v(" in the command to enable HMM filter. However, HMM filter may also filter some important structures. One should read "),e("a",{attrs:{href:"https://doi.org/10.1039/C9CP05091D",target:"_blank",rel:"noopener noreferrer"}},[a._v("the original paper"),e("OutboundLink")],1),a._v(" and carefully adjust HMM parameters.")]),a._v(" "),e("p",[a._v("HMM filter should not be used to generate the special result by time step.")]),a._v(" "),e("h3",{attrs:{id:"gui-version"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#gui-version"}},[a._v("#")]),a._v(" GUI version")]),a._v(" "),e("p",[a._v("You can open a GUI version for ReacNetGenerator by typing:")]),a._v(" "),e("div",{staticClass:"language-bash extra-class"},[e("pre",{pre:!0,attrs:{class:"language-bash"}},[e("code",[a._v("reacnetgeneratorgui\n")])])]),e("h3",{attrs:{id:"python-interface"}},[e("a",{staticClass:"header-anchor",attrs:{href:"#python-interface"}},[a._v("#")]),a._v(" Python interface")]),a._v(" "),e("p",[a._v("You can use the Python interface:")]),a._v(" "),e("div",{staticClass:"language-python extra-class"},[e("pre",{pre:!0,attrs:{class:"language-python"}},[e("code",[e("span",{pre:!0,attrs:{class:"token keyword"}},[a._v("from")]),a._v(" reacnetgenerator "),e("span",{pre:!0,attrs:{class:"token keyword"}},[a._v("import")]),a._v(" ReacNetGenerator\nReacNetGenerator"),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v("(")]),a._v("inputfiletype"),e("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),e("span",{pre:!0,attrs:{class:"token string"}},[a._v('"dump"')]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(",")]),a._v(" inputfilename"),e("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),e("span",{pre:!0,attrs:{class:"token string"}},[a._v('"dump.ch4"')]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(",")]),a._v(" atomname"),e("span",{pre:!0,attrs:{class:"token operator"}},[a._v("=")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v("[")]),e("span",{pre:!0,attrs:{class:"token string"}},[a._v("'C'")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(",")]),a._v(" "),e("span",{pre:!0,attrs:{class:"token string"}},[a._v("'H'")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(",")]),a._v(" "),e("span",{pre:!0,attrs:{class:"token string"}},[a._v("'O'")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v("]")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(")")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(".")]),a._v("runanddraw"),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v("(")]),e("span",{pre:!0,attrs:{class:"token punctuation"}},[a._v(")")]),a._v("\n")])])]),e("p",[a._v("See "),e("a",{attrs:{href:"/api/",target:"_blank"}},[a._v("Python API")]),a._v(" for details.")])])}),[],!1,null,null,null);t.default=s.exports}}]);