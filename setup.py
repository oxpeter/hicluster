#!/usr/bin/env python

from distutils.core import setup

setup(name='hicluster',
      version='2.0',
      description='hicluster package for RNA Sequencing analysis',
      author='Peter Oxley',
      author_email='oxpeter+git@gmail.com',
      url='https://github.com/oxpeter/hicluster',
      py_modules=['hicluster', 'common_path'],
      requires=['argparse','scipy','numpy','matplotlib', 'progressbar',
                'operator', 'pylab', 'mpl_toolkits', 'statsmodels',
                'multiprocessing', 'itertools', 'genomepy'
                ]
     )