const StringReplacePlugin = require("string-replace-webpack-plugin");
const TerserPlugin = require('terser-webpack-plugin');
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const HtmlWebpackPlugin = require('html-webpack-plugin');
const HtmlWebpackInlineSourcePlugin = require('html-webpack-inline-source-plugin');
const webpack = require('webpack');
const WebpackCdnPlugin = require('webpack-cdn-plugin');


const banner = `ReacNetGenerator (https://reacnetgenerator.njzjz.win/)
Copyright 2018-2019 East China Normal University
Date: ${new Date().toLocaleString()}`;

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
          'css-loader',
        {
          loader: StringReplacePlugin.replace({
            replacements: [{
              pattern: /..\/img\/bg-masthead.jpg/g,
              replacement: function (_match, _p1, _offset, _string) {
                return '';
              }
            }]
          })
        },
        {
          loader: 'postcss-loader',
          options: {
            plugins: function () {
              return [
                require('postcss-import'),
                require('precss'),
                require('cssnano')({
                  preset: 'advanced',
                }),
              ];
            }
          }
        }, {
          loader: 'sass-loader'
        }],
      },
      {
        test: /\.(jpg|png)$/,
        loader: 'url-loader',
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
      inject: 'true',
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
        removeOptionalTags: true,
        removeRedundantAttributes: true,
        removeEmptyAttributes: true,
        sortAttributes: true,
        sortClassName: true,
        useShortDoctype: true,
        ignoreCustomFragments: [/<%[\s\S]*?%>/, /<\?[\s\S]*?\?>/, /{{[\s\S]*?}}/],
        processScripts: ['text/x-jsrender']
      }
    }),
    process.env.CDN == 'yes' ? new WebpackCdnPlugin({
      modules: [
        { name: "jquery", var: ["$", "jQuery"], path: "dist/jquery.min.js" },
        { name: "regenerator-runtime", path: "runtime.min.js" },
        {
          name: "bootstrap",
          path: "dist/js/bootstrap.min.js",
          //style: "dist/css/bootstrap.min.css"
        },
        { name: "jquery.easing", path: "jquery.easing.min.js" },
        { name: "jsrender", path: "jsrender.min.js" },
        { name: "paginationjs", path: "dist/pagination.min.js" },
        {
          name: "magnific-popup",
          path: "dist/jquery.magnific-popup.min.js",
          //style: "dist/magnific-popup.css"
        },
        { nane: "bootstrap-select", path: "dist/js/bootstrap-select.min.js" },
        {
          name: "startbootstrap-creative",
          path: "js/creative.min.js",
          //style: "css/creative.min.css"
        },
        { name: "d3", path: "d3.min.js"},
        { name: "query-string", path: "index.min.js"},
        { name: "@njzjz/jsnetworkx", path: "jsnetworkx.js"},
      ],
      publicPath: "/node_modules",
      prodUrl: "//cdn.jsdelivr.net/npm/:name@:version/:path"
    }) : new HtmlWebpackInlineSourcePlugin(),
    new webpack.BannerPlugin(banner),
  ],
  optimization: {
    minimizer: [
      new TerserPlugin(),
    ],
  },
  performance: {
    hints: false,
  },
  stats: {
    entrypoints: false,
    children: false
  },
}
