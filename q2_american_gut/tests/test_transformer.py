# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import numpy as np
import qiime2

from pandas.util.testing import assert_frame_equal
from q2_american_gut.plugin_setup import QiitaMetadataFormat, QiitaMetadata
from qiime2.plugin.testing import TestPluginBase

class TestTransformers(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_pd_dataframe_to_qiita_metadata_format(self):
        transformer = self.get_transformer(pd.DataFrame, QiitaMetadataFormat)

        exp_index = pd.Index(['Sample1', 'Sample2'], dtype=object)
        exp_cols = ['Value']
        exp = pd.DataFrame([0.55, 0.99], index=exp_index, columns=exp_cols)

        obs = transformer(exp)
        obs = pd.DataFrame.from_csv(str(obs), sep='\t')
        assert_frame_equal(exp, obs)
    
    def test_qiita_metadata_format_to_pd_dataframe(self):
        filename = 'qiita-metadata.tsv'
        _, obs = self.transform_format(QiitaMetadataFormat, pd.DataFrame,
                                       filename)
        
        c1 = ['123.456.789', 'x', 'y', '1.0']
        c2 = ['123.456.012', 'thing', '0', '2.3']
        c3 = ['123.xxx.789', np.nan, 'stuff', 'None']
        c4 = ['321.xxx.789', np.nan, 'stuff', 'None']
        cols = ['#SampleID', 'foo', 'bar', 'baz']
        exp = pd.DataFrame([c1, c2, c3, c4], columns=cols)

        assert_frame_equal(exp, obs)
    
    def test_qiita_metadata_format_to_metadata(self):
        filename = 'qiita-metadata.tsv'
        _, obs = self.transform_format(QiitaMetadataFormat, qiime2.Metadata,
                                       filename)

        c1 = ['123.456.789', 'x', 'y', '1.0']
        c2 = ['123.456.012', 'thing', '0', '2.3']
        c3 = ['123.xxx.789', np.nan, 'stuff', 'None']
        c4 = ['321.xxx.789', np.nan, 'stuff', 'None']
        cols = ['#SampleID', 'foo', 'bar', 'baz']
        exp = pd.DataFrame([c1, c2, c3, c4], columns=cols)
        exp.set_index('#SampleID', inplace=True)

        exp_md = qiime2.Metadata(exp)

        self.assertEqual(obs, exp_md)
    
    def test_metadata_to_qiita_metadata_format(self):
        transformer = self.get_transformer(qiime2.Metadata, QiitaMetadataFormat)
        
        c1 = ['123.456.789', 'x', 'y', '1.0']
        c2 = ['123.456.012', 'thing', '0', '2.3']
        c3 = ['123.xxx.789', np.nan, 'stuff', 'None']
        c4 = ['321.xxx.789', np.nan, 'stuff', 'None']
        cols = ['#SampleID', 'foo', 'bar', 'baz']
        exp = pd.DataFrame([c1, c2, c3, c4], columns=cols)
        exp.set_index('#SampleID', inplace=True)

        exp_md = qiime2.Metadata(exp)
        obs = transformer(exp_md)
        obs = pd.DataFrame.from_csv(str(obs), sep='\t')
        assert_frame_equal(exp, obs)


if __name__ == '__main__':
    unittest.main()
