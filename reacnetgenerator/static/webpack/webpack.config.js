const StringReplacePlugin = require("string-replace-webpack-plugin");
const TerserPlugin = require('terser-webpack-plugin');
const MiniCssExtractPlugin = require("mini-css-extract-plugin");
const HtmlWebpackPlugin = require('html-webpack-plugin');
const HtmlWebpackInlineSourcePlugin = require('html-webpack-inline-source-plugin');

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
        test: /\.css$/,
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
				removeScriptTypeAttributes: true,
        removeStyleLinkTypeAttributes: true,
			}
    }),
    new HtmlWebpackInlineSourcePlugin()
  ],
  optimization: {
    minimizer: [
      new TerserPlugin({
        cache: false,
        parallel: true,
        terserOptions: {
          comments: false
        }
      }),
    ],
  }
}
