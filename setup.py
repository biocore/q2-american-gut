# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.develop import develop

import versioneer


def _post(obj):
    import urllib.request
    import shutil
    import os

    assets_dir = os.path.join(obj.install_libbase,
                              'q2_american_gut/assets/')

    if not os.path.exists(assets_dir):
        os.mkdir(assets_dir)

    version = versioneer.get_version()
    if 'dirty' in version:
        # we're operating in dev
        version = '2018.2'

    # source the taxonomy classifier (for assignments for deblur sequences)
    src_url = 'https://data.qiime2.org/%s/common/gg-13-8-99-515-806-nb-classifier.qza' % version  # noqa
    out_f = os.path.join(assets_dir, 'gg-13-8-99-515-806-nb-classifier.qza')

    req = urllib.request.Request(str(src_url), headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as response, open(out_f, 'wb') as out:
        shutil.copyfileobj(response, out)

    # source the Greengenes tree (for closed reference data)
    src_url = 'ftp://ftp.microbio.me/greengenes_release/gg_13_8_otus/trees/97_otus.tree'  # noqa
    out_f = os.path.join(assets_dir, '97_otus.tree')

    req = urllib.request.Request(str(src_url), headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as response, open(out_f, 'wb') as out:
        shutil.copyfileobj(response, out)


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        install.run(self)
        _post(self)


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)
        _post(self)


cmdclass = versioneer.get_cmdclass()
cmdclass['install'] = PostInstallCommand
cmdclass['develop'] = PostDevelopCommand

setup(
    name="q2-american-gut",
    version=versioneer.get_version(),
    packages=find_packages(),
    package_data={'q2_american_gut.tests': ['data/*'],
                  'q2_american_gut': ['assets/reference_phylogeny_tiny.qza',
                                      'assets/reference_alignment_tiny.qza',
                                      'assets/report/index.html',
                                      'assets/report/resources/*']
                  },
    author="Daniel McDonald",
    author_email="danielmcdonald@ucsd.edu",
    description="American Gut processing and interaction",
    license='BSD-3-Clause',
    url="http://americangut.org",
    cmdclass=cmdclass,
    entry_points={
        'qiime2.plugins': ['q2-american-gut=q2_american_gut.plugin_setup:'
                           'plugin']
    },
    zip_safe=False,
)
