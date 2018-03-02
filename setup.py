#!/usr/bin/env python3
"""
Copyright 2018 Jane Hawkey (jane.hawkey@unimelb.edu.au)
https://github.com/jhawkey/sra_read_downloader

Setup script for SRA Read Downloader

This file is part of SRA Read Downloader. SRA Read Downloader is free software: you can
redistribute it and/or modify it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License, or (at your option) any later
version. SRA Read Downloader is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU General Public License for more details. You should have received a copy of
the GNU General Public License along with SRA Read Downloader. If not, see
<http://www.gnu.org/licenses/>.
"""

from setuptools import setup

with open('README.md', 'rb') as readme:
    long_description = readme.read()
if not isinstance(long_description, str):
    long_description = long_description.decode()

# Get the version from sra_read_downloader.
__version__ = '0.0.0'
script_lines = open('sra_read_downloader.py').readlines()
version_line = [x for x in script_lines if x.startswith('__version__')][0]
exec(version_line)

setup(name='SRA Read Downloader',
      version=__version__,
      description='A handy script for downloading reads from any kind of accession from the SRA',
      long_description=long_description,
      url='https://github.com/jhawkey/sra_read_downloader',
      author='Jane Hawkey',
      author_email='jane.hawkey@unimelb.edu.au',
      license='GPLv3',
      scripts=['sra_read_downloader.py'],
      install_requires=['pandas'],
      zip_safe=False)
