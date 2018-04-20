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


if __name__ == '__main__':
    unittest.main()

