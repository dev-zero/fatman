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

# Create a secret
python -c "import os; import binascii; print('SECRET_KEY = {}'.format(binascii.hexlify(os.urandom(24))));" > fatman.cfg

# Run the (development) server
FATMAN_SETTINGS=$PWD/fatman.cfg ./manage.py runserver

# Open the interface
xdg-open http://localhost:5000/admin/
```
