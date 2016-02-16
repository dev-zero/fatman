
from setuptools import setup

setup(
    name='fatman',
    version='0.1.dev0',
    packages=['fatman'],
    license='GPL3',
    install_requires=[
        "Flask",
        "peewee",
        "psycopg2",
        "Flask-RESTful",
        "Flask-Admin",
        "wtf-peewee",
        "Flask-Script",
        ]
    )
