# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from qiime2.plugin import Plugin, Metadata, List, Str
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.ordination import PCoAResults
from ._visualizer import report

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


plugin.visualizers.register_function(
    function=report,
    inputs={
        'taxonomy': FeatureData[Taxonomy],
        'table': FeatureTable[Frequency],
        'pcoa': PCoAResults,
        'alpha': SampleData[AlphaDiversity]
        },
    parameters={'metadata': Metadata,
                'samples': List[Str]},
    input_descriptions={
        'taxonomy': 'placeholder 1',
        'table': 'Feature table to visualize at various taxonomic levels.',
        'pcoa': 'placeholder',
        'alpha': 'placeholder'},
    parameter_descriptions={'metadata': 'The sample metadata.',
                            'samples': 'list of relevant samples.'},
    name='Generate AGP report.',
    description='none'
)
