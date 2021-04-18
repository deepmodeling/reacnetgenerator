const StringReplacePlugin = require("string-replace-webpack-plugin");
const TerserPlugin = require('terser-webpack-plugin');
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const HtmlWebpackPlugin = require('html-webpack-plugin');
const ScriptExtHtmlWebpackPlugin = require("script-ext-html-webpack-plugin");
const HTMLInlineCSSWebpackPlugin = require("html-inline-css-webpack-plugin").default;
const OptimizeCssAssetsPlugin = require('optimize-css-assets-webpack-plugin');
const webpack = require('webpack');

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
                  require('precss'),
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
        ignoreCustomFragments: [ /<%[\s\S]*?%>/, /<\?[\s\S]*?\?>/, /{{[\s\S]*?}}/ ],
        processScripts: ['text/x-jsrender']
			}
    }),
	new HTMLInlineCSSWebpackPlugin(),
	new ScriptExtHtmlWebpackPlugin({inline: /\.js$/}),
    new webpack.BannerPlugin(banner)
  ],
  optimization: {
    minimizer: [
      new TerserPlugin(),
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
}
