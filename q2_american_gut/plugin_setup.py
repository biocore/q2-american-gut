# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from qiime2.plugin import Plugin
from q2_types.feature_table import FeatureTable, Frequency

import q2_american_gut


plugin = Plugin(
    name='american-gut',
    version=q2_american_gut.__version__,
    website='https://github.com/biocore/q2-american-gut',
    package='q2_american_gut',
    description=('This QIIME 2 plugin supports processing and utilizing '
                 'American Gut Project data'),
    short_description='Plugin for exploring American Gut data.',
    citation_text='https://doi.org/10.1101/277970'
)

def dummy(foo: biom.Table) -> biom.Table:
    return foo

plugin.methods.register_function(
    function=dummy,
    inputs={'foo': FeatureTable[Frequency]},
    parameters={},
    outputs=[('bar', FeatureTable[Frequency]), ],
    input_descriptions={
        'foo': "Same same"
    },
    output_descriptions={'bar': 'Same same'},
    name='A dummy placeholder',
    description="A dummy placeholder"
)
