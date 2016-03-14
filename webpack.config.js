var webpack = require('webpack');

module.exports = {
    entry: {
        'app': './fatman/app.ts',
        'vendor': './fatman/vendor.ts'
    },
    output: {
        path: __dirname + '/fatman/static',
        filename: 'bundle.js'
    },
    resolve: {
        extensions: ['', '.js', '.ts']
    },
    devtool: 'source-map',
    module: {
        loaders: [
            {
                test: /\.ts/,
                loaders: ['ts-loader'],
                exclude: /node_modules/
            }
        ]
    },
    plugins: [
        new webpack.optimize.CommonsChunkPlugin('vendor', 'vendor.bundle.js')
    ]
}
