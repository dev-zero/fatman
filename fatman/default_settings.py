
SQLALCHEMY_DATABASE_URI = 'postgresql:///fatman'
SQLALCHEMY_TRACK_MODIFICATIONS = False

APPLICATION_ROOT = ''
DEBUG = False

# File-Upload defaults
UPLOADED_RESULTS_DEST = '/var/empty'
UPLOADED_RESULTS_ALLOW = ['md', 'json', 'xml', 'zip', 'tar', 'tgz', 'gz', 'tbz2', 'bz2', 'xz']

# Flask-Security defaults
SECURITY_PASSWORD_HASH = 'pbkdf2_sha512'
SECURITY_PASSWORD_SALT = '' # don't worry, read https://github.com/mattupstate/flask-security/issues/470
