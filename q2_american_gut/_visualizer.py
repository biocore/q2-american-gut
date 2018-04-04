# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources
import shutil
import skbio
import biom

import pandas as pd
import q2templates

from qiime2 import Metadata


TEMPLATES = pkg_resources.resource_filename('q2_american_gut', 'assets')


def report(output_dir: str,
           pcoa: skbio.OrdinationResults,
           metadata: Metadata,
           alpha: pd.Series,
           table: biom.Table,
           taxonomy: pd.Series,
           samples: list) -> None:
    metadata = metadata.to_dataframe()

    _insanity_checker(samples, metadata, table, alpha, pcoa)

    index = os.path.join(TEMPLATES, 'report', 'index.html')
    q2templates.render(index, output_dir, context={'name': 'foo'})

    # Copy assets for rendering figure
    shutil.copytree(os.path.join(TEMPLATES, 'report', 'resources'),
                    os.path.join(output_dir, 'resources'))


def _insanity_checker(samples, metadata, table, alpha, pcoa):
    samples = set(samples)

    if not samples.issubset(set(metadata.index)):
        raise ValueError('There are missing samples in the metadata')
    if not samples.issubset(set(table.ids('sample'))):
        raise ValueError('There are missing samples in the BIOM table')
    if not samples.issubset(alpha.index):
        raise ValueError('There are missing samples in the alpha diversity '
                         'vector')
    if not samples.issubset(set(pcoa.samples.index)):
        raise ValueError('There are missing samples in the ordination')
