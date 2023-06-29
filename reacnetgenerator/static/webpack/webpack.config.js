// SPDX-License-Identifier: LGPL-3.0-or-later
const StringReplacePlugin = require("string-replace-webpack-plugin");
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const HtmlWebpackPlugin = require('html-webpack-plugin');
const ScriptExtHtmlWebpackPlugin = require("script-ext-html-webpack-plugin");
const HTMLInlineCSSWebpackPlugin =
    require("html-inline-css-webpack-plugin").default;
const OptimizeCssAssetsPlugin = require('css-minimizer-webpack-plugin');
const WebpackCdnPlugin = require('webpack-cdn-plugin');
const webpack = require('webpack');

const year = new Date().getFullYear();
const banner = `ReacNetGenerator (https://reacnetgenerator.njzjz.win/)
Copyright 2018-${year} East China Normal University
Date: ${new Date().toLocaleString()}`;
const buildweb = process.env.REACNETGENERATOR_BUILDWEB;

module.exports = {
  entry: __dirname + "/reacnetgen.js",
  output: {
    path: __dirname + "/",
    filename: "bundle.js"
  },
  mode: 'production',
  module: {
    rules: [
      {
        test: /\.scss$/,
        use: [{
          loader: MiniCssExtractPlugin.loader,
        },
        {
          loader:'css-loader',
          options:{
            url: !buildweb,
          }
        },
        {
          loader: StringReplacePlugin.replace({
            replacements: [{
              pattern: /..\/assets\/img\/bg-masthead.jpg/g,
              replacement: function (_match, _p1, _offset, _string) {
                return '';
              }
            }]
          })
        },
        {
          loader: 'postcss-loader',
          options: {
            postcssOptions: {
              plugins: function () {
                return [
                  require('postcss-import'),
                  require('cssnano')({
                    preset: 'advanced',
                  }),
                ];
              }
            }
          }
        }, {
          loader: 'sass-loader'
        }],
      },
      {
        test: /\.(jpg|png)$/,
        type: 'asset/inline',
      },
      {
        test: /\.js$/,
        exclude: /node_modules/,
        use: [
          { loader: "ifdef-loader" }
        ]
      },
    ],
  },
  plugins: [
    new StringReplacePlugin(),
    new MiniCssExtractPlugin({
      filename: "bundle.css"
    }),
    new HtmlWebpackPlugin({
      filename: 'bundle.html',
      template: __dirname + '/template.html',
      inject: 'body',
      inlineSource: '.(js|css)$',
      minify: {
        collapseWhitespace: true,
        collapseBooleanAttributes: true,
        collapseInlineTagWhitespace: true,
        minifyCSS: true,
        minifyJS: true,
        minifyURLs: true,
        decodeEntities: true,
        removeScriptTypeAttributes: true,
        removeStyleLinkTypeAttributes: true,
        removeAttributeQuotes: true,
        removeOptionalTags: false,
        removeRedundantAttributes: true,
        removeEmptyAttributes: true,
        sortAttributes: true,
        sortClassName: true,
        useShortDoctype: true,
        ignoreCustomFragments: [/<%[\s\S]*?%>/, /<\?[\s\S]*?\?>/, /{{[\s\S]*?}}/],
        processScripts: ['text/x-jsrender']
      }
    }),
    ...(buildweb ? [
      // build web, replace CDN
      new WebpackCdnPlugin({
        modules: [
          { name: "jquery", var: "$", path: "dist/jquery.min.js" },
          { name: "regenerator-runtime", path: "runtime.js", var: "regeneratorRuntime" },
          { name: "@popperjs/core", path: "dist/umd/popper.min.js", var: "Popper" }, // must above bootstrap
          {
            name: "bootstrap",
            path: "dist/js/bootstrap.min.js",
            style: "dist/css/bootstrap.min.css"
          },
          { name: "jsrender", path: "jsrender.min.js", var: "$.jsrender" },
          {
            name: "paginationjs",
            path: "dist/pagination.min.js",
            style: "dist/pagination.css",
            var: "$.paginationjs"
          },
          {
            name: "magnific-popup",
            path: "dist/jquery.magnific-popup.min.js",
            style: "dist/magnific-popup.css",
            var: "$.magnificPopup"
          },
          {
            name: "bootstrap-select",
            path: "dist/js/bootstrap-select.min.js",
            style: "dist/css/bootstrap-select.min.css",
            var: "$.selectpicker"
          },
          {
            name: "startbootstrap-creative",
            path: "dist/js/scripts.js",
            style: "dist/css/styles.css"
          },
          { name: "d3", path: "d3.min.js" },
          { name: "@njzjz/jsnetworkx", path: "jsnetworkx.js", var: "jsnx" },
          { name: "animejs", path: "lib/anime.min.js", var: "anime" },
          { name: "smiles-drawer", path: "dist/smiles-drawer.min.js", var: "SmilesDrawer"},
        ],
        prodUrl: '//unpkg.com/:name@:version/:path',
      })
    ] : [
      // build inline
      new HTMLInlineCSSWebpackPlugin(),
      new ScriptExtHtmlWebpackPlugin({ inline: /\.js$/ }),
    ]),
    new webpack.BannerPlugin(banner)
  ],
  optimization: {
    minimizer: [
      `...`,
      new OptimizeCssAssetsPlugin(),
    ],
  },
  performance: {
    hints: false,
  },
  stats: {
    entrypoints: false,
    children: false
  },
};
