# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib
from qiime2.plugin import Plugin, Int, Str, Choices, Metadata, Bool, List
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy
from q2_types.sample_data import AlphaDiversity, SampleData
from q2_types.ordination import PCoAResults
from q2_types.tree import Phylogeny, Rooted

import q2_american_gut
from q2_american_gut._visualizer import report
from q2_american_gut._type import QiitaMetadata
from q2_american_gut._format import (QiitaMetadataFormat,
                                     QiitaMetadataDirectoryFormat)


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

plugin.register_formats(QiitaMetadataFormat, QiitaMetadataDirectoryFormat)

plugin.register_semantic_types(QiitaMetadata)
plugin.register_semantic_type_to_format(
    QiitaMetadata,
    artifact_format=QiitaMetadataDirectoryFormat
)


# TODO: add support for shotgun retrieval
# TODO: add support for metabolomic retrieval
# TODO: add support for HMP reference genome hits
plugin.methods.register_function(
    function=q2_american_gut.fetch_amplicon,
    name='Fetch amplicon data',
    description=('This method obtains study amplicon data from Qiita.'),
    inputs={},
    input_descriptions={},
    parameters={
        'qiita_study_id': Str,
        'processing_type': Str % Choices(['deblur', 'closed-reference']),
        'trim_length': Str % Choices(['90', '100', '150']),
        'threads': Int,
        'debug': Bool
    },
    parameter_descriptions={
        'qiita_study_id': 'The study to obtain',
        'processing_type': 'How the OTUs were assessed',
        'trim_length': 'The sequence trim length to use',
        'threads': ('Number of parallel downloads to perform.'),
        'debug': ('Whether to operate in debug mode. If debug mode, a small '
                  'subset of data are fetched.')
    },
    outputs=[
        ('feature_table', FeatureTable[Frequency]),
        ('feature_taxonomy', FeatureData[Taxonomy]),
        ('sample_metadata', QiitaMetadata),
        ('phylogeny', Phylogeny[Rooted])
    ],
    output_descriptions={
        'feature_table': "A feature table of the sample data",
        'feature_taxonomy': "Feature taxonomy information",
        'sample_metadata': "Feature metadata",
        'phylogeny': "A phylogeny relating the features"
    }
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
importlib.import_module('q2_american_gut._transformer')
