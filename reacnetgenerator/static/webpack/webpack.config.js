var StringReplacePlugin = require("string-replace-webpack-plugin");

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
        use: ['style-loader', 'css-loader',
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
    new StringReplacePlugin()
  ]
}
