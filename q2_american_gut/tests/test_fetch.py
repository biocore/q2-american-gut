# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import biom
import numpy as np
import skbio
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAIterator

from q2_american_gut import fetch_amplicon
from q2_american_gut._fetch import _determine_context, _get_featuredata_from_table



class DNAIteratorTests(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_empty_table(self):
        tab = biom.Table([], [], [])

        with self.assertRaisesRegex(ValueError, "No features"):
            _get_featuredata_from_table(tab)

    def test_get_iterator(self):
        tab = biom.Table(np.ones((3, 2)), ['ATCC', 'ATGG', 'CACA'],
                         ['S1', 'S2'])
        exp = [skbio.DNA('ATCC', metadata={'id': 'ATCC'}),
               skbio.DNA('ATGG', metadata={'id': 'ATGG'}),
               skbio.DNA('CACA', metadata={'id': 'CACA'})]
        obs = list(_get_featuredata_from_table(tab))
        self.assertEqual(obs, exp)


class DetermineContextTests(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_context_without_processing_type(self):
        proctype = 'something wrong'
        trim = 113

        msg = "Cannot find %s-%dnt for 16S data" % (proctype, trim)
        with self.assertRaisesRegex(ValueError, msg):
            _determine_context(proctype, trim)

    def test_context_does_not_have_trim_length(self):
        proctype = 'deblur'
        trim = 42

        msg = "Cannot find %s-%dnt for 16S data" % (proctype, trim)
        with self.assertRaisesRegex(ValueError, msg):
            _determine_context(proctype, trim)

    def test_valid_deblur_context(self):
        proctype = 'deblur'
        trim = 100
        ctx = _determine_context(proctype, trim)
        self.assertIn('deblur', ctx.lower())
        self.assertIn('100nt', ctx.lower())

    def test_valid_closed_reference_context(self):
        proctype = 'closed-reference'
        trim = 100
        ctx = _determine_context(proctype, trim)
        self.assertIn('closed-reference', ctx.lower())
        self.assertIn('100nt', ctx.lower())


class TestFetch(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_non_existant_study_id(self):
        id_ = "99999999999999"
        with self.assertRaisesRegex(ValueError, "has no samples"):
            fetch_amplicon(id_, 'deblur', 100)

    def test_study_id_must_be_positive_integer(self):
        id_ = "-1"
        with self.assertRaisesRegex(ValueError,
                                    "ID %s is not a qiita ID" % id_):
            fetch_amplicon(id_, 'deblur', 100)

        id_ = "assd"
        with self.assertRaisesRegex(ValueError,
                                    "ID %s is not a qiita ID" % id_):
            fetch_amplicon(id_, 'deblur', 100)


    def test_positive_fetch_deblur(self):
        obs_table, obs_tax, obs_md, obs_phy = \
                            fetch_amplicon('2136', 'deblur', 100)
    
        view_table = obs_table.view(biom.Table)
        view_tax = obs_tax.view(pd.DataFrame)
        view_md = obs_md.view(pd.DataFrame)
        exp_table_shape = (15448, 504)
        self.assertEqual(view_table.shape, exp_table_shape)
        # view taxonomy as pd.Dataframe, assert that all identifiers are present
        # for md assert that all 
        # for phyl want to test that the feature identifiers in table are subest of 
        # the tip names in obs_phy are is subset of tips in table 
        df_table = view_table.to_dataframe()
        self.assertEqual(set(df_table.index), set(view_md.index))
        

    def test_positive_fetch_closed_ref(self):
        obs_table, obs_tax, obs_md, obs_phy = fetch_amplicon('2136', 
                                            'closed-reference',  100)

        view_table = obs_table.view(biom.Table)
        view_tax = obs_tax.view(pd.DataFrame)
        view_md = obs_md.view(pd.DataFrame)
        exp_table_shape = (7432, 533)
        self.assertEqual(view_table.shape, exp_table_shape)
        df_table = view_table.to_dataframe()
        self.assertEqual(set(df_table.index), set(view_md.index))
        
    def test_get_closed_reference_study(self):
        id_ = '10343'
        proc_type = 'closed-reference'
        length = 100
        debug = True

        table, tax, md, tree = fetch_amplicon(id_, proc_type,
                                              length, debug=debug)

        exp_ids = ['10343.1384a.36263',
                   '10343.2024a.36263',
                   '10343.2025a.36263',
                   '10343.2026a.36263',
                   '10343.2160a.36263',
                   '10343.2161a.36263',
                   '10343.2162a.36263',
                   '10343.BLANK.JS6.12E.36263',
                   '10343.BLANK.JS6.12F.36263',
                   '10343.BLANK.JS6.12G.36263']

        # test consistency between the outputs
        self.assertEqual(sorted(table.ids()), sorted(exp_ids))
        self.assertEqual(len(table.ids(axis='observation')), 826)
        self.assertEqual(set(table.ids(axis='observation')),
                         set(tax.index))
        self.assertTrue(set(tax.index).issubset({n.name for n in tree.tips()}))
        self.assertEqual(set(table.ids()), set(md.index))

        
        def test_fetch_phylogeny(self):
            
            pass

        def test_fetch_taxonomy(self):

            pass

if __name__ == '__main__':
    unittest.main()

