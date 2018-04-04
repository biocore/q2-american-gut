# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import versioneer

setup(
    name="q2-american-gut",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'q2_american_gut.tests': ['data/*'],
                  'q2_american_gut/': ['assets/report/index.html',
                                       'assets/report/resources/*']
                  },
    author="Daniel McDonald",
    author_email="danielmcdonald@ucsd.edu",
    description="American Gut processing and interaction",
    license='BSD-3-Clause',
    url="http://americangut.org",
    entry_points={
        'qiime2.plugins': ['q2-american-gut=q2_american_gut.plugin_setup:'
                           'plugin']
    },
    zip_safe=False,
)
