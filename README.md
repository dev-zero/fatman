# FATMAN - Flexible Accuracy Testing MANager

## How to run

First create a (PostgreSQL) database named `fatman` and make sure the user running the server has access to it.

Then:

```sh
# Grab the code:
git clone ...

# Change the directory
cd fatman/

# Create a virtual environment and use it
virtualenv venv
source venv/bin/activate

# Install required packages
pip install -r requirements.txt

# Create an initial configuration (including a custom secret key for the session cookies, will overwrite an existing `fatman.cfg` in your current directory)
./manage.py createconfig

# Initialize the database (this is idempotent)
./manage.py initdb

# Run the (development) server
FATMAN_SETTINGS=$PWD/fatman.cfg ./manage.py runserver

# Open the interface
xdg-open http://localhost:5000/admin/
```

## Web-Application Deployment

The web-frontend is based on AngularJS 2.0 + Bootstrap 3 and written in TypeScript.

To build/deploy it you need the following packages on your system:

* NodeJS
* NPM

After that, you can bootstrap it as follows:

```sh
# let NPM install all required JavaScript/TypeScript packages (in ./node_modules)
npm install

# setup required typings for TypeScript
./node_modules/.bin/typings install

# Build the JavaScript bundles
./node_modules/.bin/webpack
```

After that, the web-application should be available on:

  http://localhost:5000/app
