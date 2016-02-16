
from setuptools import setup

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='fatman',
    version='0.1.dev0',
    packages=['fatman'],
    license='GPL3',
    install_requires=requirements,
    )
