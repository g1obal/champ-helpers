"""
CHAMP Helpers setup

Author: Gokhan Oztarhan
champ-helpers created date: 09/06/2019
Setup module created date: 28/12/2021
Setup module last modified: 26/11/2022
"""

from distutils.core import setup
import re

from champio import __version__


# Parse long_description.
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Parse requirements. Keep only package names. The following are the minimal
# requirements and do not include version numbers (see README.md).
with open('requirements.txt', 'r', encoding='utf-8') as f:
    delimeters = '<|<=|==|!=|~=|>=|>'
    requirements = [
        re.split(delimeters, line)[0].strip() \
        for line in f.readlines() if not line.startswith('#')
    ]

# Initialize setup
setup(
    name = 'champ-helpers',
    version = __version__,
    description = 'A set of helper modules and scripts for CHAMP.',
    long_description = long_description,
    author = 'Gokhan Oztarhan',
    author_email = 'gooztarhan@gmail.com',
    url = 'https://github.com/g1obal/champ-helpers',
    license = 'MIT License',
    license_files = ['LISENCE'],
    keywords = [
        'python3', 'CHAMP', 'mpi', 'data', 'data analysis', 'plotting',
    ],
    packages = [
        'champio', 
        'champio.latgen', 
    ],
    package_dir = {
        'champio': 'champio', 
        'champio.latgen': 'champio/latgen',
    },
    install_requires = requirements,
)

