# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd

import qiime2
from pandas.util.testing import assert_series_equal
from q2_american_gut.plugin_setup import QiitaMetadataFormat, QiitaMetadata
from qiime2.plugin.testing import TestPluginBase

# want metdata to format, format to metadata, qiitameta to meta and meta to qiita
class TestTransformers(TestPluginBase):
    package = "q2_american_gut.tests"
    
    def test_pd_dataframe_to_qiita_metadata_format(self):
        # this transformer i believe needs to be created :L
        transformer = self.get_transformer(pd.Series, QiitaMetadataFormat)

        exp_index = 'some indexes...'
        exp = 'some dataframe pd.DataFrame(...)'

        obs = transformer(exp)
        obs = pd.Series.from_csv(str(obs), sep='\t', header=0)

        assert_series_equal(exp, obs)
        
    def test_qiita_metadata_format_to_pd_dataframe(self):
        filename = 'qiita-metadata.tsv'
        _, obs = self.transform_format(QiitaMetadataFormat, pd.DataFrame,
                                       filename)
        exp_index = 'some indexes pd.Index(...)'
        exp = 'some dataframe pd.DataFrame(...)'
        assert_series_equal(exp, obs)

    def test_qiita_metadata_format_to_metadata(self):
        filename = 'qiita-metadata.tsv'
        _, obs = self.transform_format(QiitaMetadataFormat, qiime2.Metadata,
                                       filename)

        exp_index = 'some indexes? pd.Index(...)'
        exp_df = 'some dataframe pd.DataFrame(...)'
        exp_md = qiime2.Metadata(exp_df)

        self.assertEqual(obs, exp_md)

    def test_metadata_to_qiita_metadata_format(self):
        # i believe this is is how you would test this transoformation, 
        # using a transformer as in pd_dataframe_to_qiita_md_format, 
        # but im not sure. If this is the case, then this transformer also
        # needs to be created
        transformer = self.get_transformer(qiime2.Metadata, QiitaMetadataFormat)
        # this is unfinished
        pass

    def test_metadata_to_qiita_metadata(self):
        pass # this isnt even started

    def test_qiita_metadata_to_metadata(self):
        pass


if __name__ == '__main__':
    unittest.main()
