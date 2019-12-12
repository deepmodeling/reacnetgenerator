(window.webpackJsonp=window.webpackJsonp||[]).push([[5],{205:function(a,e,t){"use strict";t.r(e);var n=t(0),r=Object(n.a)({},(function(){var a=this,e=a.$createElement,t=a._self._c||e;return t("ContentSlotsDistributor",{attrs:{"slot-key":a.$parent.slotKey}},[t("video",{attrs:{width:"512",height:"397.6",controls:""}},[t("source",{attrs:{src:"http://www.rsc.org/suppdata/c9/cp/c9cp05091d/c9cp05091d2.mp4",type:"video/mp4"}})]),a._v(" "),t("h1",{attrs:{id:"citation-and-contact"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#citation-and-contact"}},[a._v("#")]),a._v(" Citation and contact")]),a._v(" "),t("p",[a._v("ReacNetGenerator: an Automatic Reaction Network Generator for Reactive Molecular Dynamic Simulations, Phys. Chem. Chem. Phys., 2019, doi: "),t("a",{attrs:{href:"https://dx.doi.org/10.1039/C9CP05091D",target:"_blank",rel:"noopener noreferrer"}},[a._v("10.1039/C9CP05091D"),t("OutboundLink")],1)]),a._v(" "),t("p",[a._v("jinzhe.zeng@rutgers.edu (Jinzhe Zeng), tzhu@lps.ecnu.edu.cn (Tong Zhu)")]),a._v(" "),t("h1",{attrs:{id:"installation"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#installation"}},[a._v("#")]),a._v(" Installation")]),a._v(" "),t("p",[a._v("First, you need to download the source code from "),t("a",{attrs:{href:"https://github.com/tongzhugroup/reacnetgenerator/releases",target:"_blank",rel:"noopener noreferrer"}},[a._v("the Releases page"),t("OutboundLink")],1),a._v(". Then install ReacNetGenerator with one of the following guides:")]),a._v(" "),t("h2",{attrs:{id:"building-a-conda-package"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#building-a-conda-package"}},[a._v("#")]),a._v(" Building a conda package")]),a._v(" "),t("ol",[t("li",[t("a",{attrs:{href:"https://conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html",target:"_blank",rel:"noopener noreferrer"}},[a._v("Install Anaconda or Miniconda"),t("OutboundLink")],1),a._v(" to obtain conda.")]),a._v(" "),t("li",[a._v("Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:")])]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("conda config --add channels conda-forge\nconda build conda/recipe\nconda "),t("span",{pre:!0,attrs:{class:"token function"}},[a._v("install")]),a._v(" reacnetgenerator --use-local\nreacnetgenerator -h\n")])])]),t("h2",{attrs:{id:"building-a-docker-image"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#building-a-docker-image"}},[a._v("#")]),a._v(" Building a Docker Image")]),a._v(" "),t("ol",[t("li",[t("a",{attrs:{href:"https://docs.docker.com/install/",target:"_blank",rel:"noopener noreferrer"}},[a._v("Install Docker"),t("OutboundLink")],1),a._v(".")]),a._v(" "),t("li",[a._v("Decompress reacnetgenerator.zip and build in the main directory of ReacNetGenerator:")])]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("docker build "),t("span",{pre:!0,attrs:{class:"token builtin class-name"}},[a._v(".")]),a._v(" -t njzjz/reacnetgenerator\ndocker run njzjz/reacnetgenerator reacnetgenerator -h\n")])])]),t("h2",{attrs:{id:"installing-via-pip"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#installing-via-pip"}},[a._v("#")]),a._v(" Installing via pip")]),a._v(" "),t("ol",[t("li",[a._v("Install "),t("a",{attrs:{href:"https://github.com/openbabel",target:"_blank",rel:"noopener noreferrer"}},[a._v("OpenBabel"),t("OutboundLink")],1),a._v(", "),t("a",{attrs:{href:"https://github.com/rdkit/rdkit",target:"_blank",rel:"noopener noreferrer"}},[a._v("RDKit"),t("OutboundLink")],1),a._v(", and "),t("a",{attrs:{href:"https://github.com/yarnpkg/yarn",target:"_blank",rel:"noopener noreferrer"}},[a._v("Yarn"),t("OutboundLink")],1),a._v(".")]),a._v(" "),t("li",[a._v("Decompress reacnetgenerator.zip and use "),t("code",[a._v("pip")]),a._v(" to install in the main directory of ReacNetGenerator. Note that a C/C++ compiler must be installed.")])]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("pip "),t("span",{pre:!0,attrs:{class:"token function"}},[a._v("install")]),a._v(" "),t("span",{pre:!0,attrs:{class:"token builtin class-name"}},[a._v(".")]),a._v("\n")])])]),t("h1",{attrs:{id:"usage"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#usage"}},[a._v("#")]),a._v(" Usage")]),a._v(" "),t("h2",{attrs:{id:"command-line"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#command-line"}},[a._v("#")]),a._v(" Command line")]),a._v(" "),t("p",[a._v("ReacNetGenerator can process any kind of trajectory files containing atomic coordinates, e.g. a LAMMPS dump file prepared by running “dump 1 all custom 100 dump.reaxc id type x y z” in LAMMPS:")]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("reacnetgenerator --dump -i dump.reaxc -a C H O\n")])])]),t("p",[a._v("where C, H, and O are atomic names in the input file. "),t("a",{attrs:{href:"/report.html?jdata=https%3A%2F%2Fgist.githubusercontent.com%2Fnjzjz%2Fe9a4b42ceb7d2c3c7ada189f38708bf3%2Fraw%2F83d01b9ab1780b0ad2d1e7f934e61fa113cb0f9f%2Fmethane.json",target:"_blank"}},[a._v("Analysis report")]),a._v(" will be generated automatically.")]),a._v(" "),t("p",[a._v("Also, ReacNetGenerator can process files containing bond information, e.g. LAMMPS bond file:")]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("reacnetgenerator -i bonds.reaxc -a C H O\n")])])]),t("p",[a._v("You can running the following script for help:")]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("reacnetgenerator -h\n")])])]),t("h2",{attrs:{id:"gui-version"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#gui-version"}},[a._v("#")]),a._v(" GUI version")]),a._v(" "),t("p",[a._v("You can open a GUI version for ReacNetGenerator by typing:")]),a._v(" "),t("div",{staticClass:"language-bash extra-class"},[t("pre",{pre:!0,attrs:{class:"language-bash"}},[t("code",[a._v("reacnetgeneratorgui\n")])])]),t("h1",{attrs:{id:"awards"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#awards"}},[a._v("#")]),a._v(" Awards")]),a._v(" "),t("ul",[t("li",[a._v("The First Prize in 2019 (the 11th Session) Shanghai Computer Application Competition for College Students")]),a._v(" "),t("li",[a._v("The First Prize in 2019 (the 12th Session) Chinese Computer Design Competition for College Students")])]),a._v(" "),t("h1",{attrs:{id:"acknowledge"}},[t("a",{staticClass:"header-anchor",attrs:{href:"#acknowledge"}},[a._v("#")]),a._v(" Acknowledge")]),a._v(" "),t("ul",[t("li",[a._v("National Natural Science Foundation of China (Grants No. 91641116)")]),a._v(" "),t("li",[a._v("National Innovation and Entrepreneurship Training Program for Undergraduate (201910269080)")]),a._v(" "),t("li",[a._v("ECNU Multifunctional Platform for Innovation (No. 001)")])])])}),[],!1,null,null,null);e.default=r.exports}}]);