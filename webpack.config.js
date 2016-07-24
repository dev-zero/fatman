var webpack = require('webpack');
var path = require('path');

// Webpack Config
var webpackConfig = {
  entry: {
    'polyfills': './fatman/app/polyfills.ts',
    'vendor': './fatman/app/vendor.ts',
    'app': './fatman/app/main.ts',
  },

  output: {
    path: path.join(__dirname, 'fatman', 'static'),
  },

  plugins: [
    new webpack.optimize.CommonsChunkPlugin({ name: ['app', 'vendor', 'polyfills'], minChunks: Infinity }),
  ],

  module: {
    loaders: [
      // .ts files for TypeScript
    { test: /\.ts/, loaders: ['ts-loader'], exclude: /node_modules/ },
    ]
  }
};


// Our Webpack Defaults
var defaultConfig = {

  devtool: 'cheap-module-source-map',
  cache: true,
  debug: true,

  output: {
    filename: '[name].bundle.js',
    sourceMapFilename: '[name].map',
    chunkFilename: '[id].chunk.js'
  },

  resolve: {
    root: [ path.join(__dirname, 'fatman') ],
    extensions: ['', '.js', '.ts'],
    alias: {
      'angular2/testing': path.join(__dirname, 'node_modules', '@angular', 'core', 'testing.js'),
      '@angular/testing': path.join(__dirname, 'node_modules', '@angular', 'core', 'testing.js'),
      'angular2/core': path.join(__dirname, 'node_modules', '@angular', 'core', 'index.js'),
      'angular2/platform/browser': path.join(__dirname, 'node_modules', '@angular', 'platform-browser', 'index.js'),
      'angular2/testing': path.join(__dirname, 'node_modules', '@angular', 'testing', 'index.js'),
      'angular2/router': path.join(__dirname, 'node_modules', '@angular', 'router-deprecated', 'index.js'),
      'angular2/http': path.join(__dirname, 'node_modules', '@angular', 'http', 'index.js'),
      'angular2/http/testing': path.join(__dirname, 'node_modules', '@angular', 'http', 'testing.js')
    },
  },

  module: {
    preLoaders: [{
      test: /\.js$/,
      loader: 'source-map-loader',
      exclude: [
        // these packages have problems with their sourcemaps
        path.join(__dirname, 'node_modules', 'rxjs'),
        path.join(__dirname, 'node_modules', '@angular2-material'),
        path.join(__dirname, 'node_modules', '@angular'),
      ]
    }],
    noParse: [
      path.join(__dirname, 'node_modules', 'zone.js', 'dist'),
      path.join(__dirname, 'node_modules', 'angular2', 'bundles')
    ]
  },
}

// Merging this Webpack config and the defaults
var webpackMerge = require('webpack-merge');
module.exports = webpackMerge(defaultConfig, webpackConfig);

// vim: set ts=2 sw=2 tw=0 :
