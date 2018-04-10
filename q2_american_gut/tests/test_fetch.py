# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_american_gut import fetch_amplicon


class TestFetch(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_non_existant_study_id(self):
        id_ = "-1"
        with self.assertRaisesRegex(ValueError, "study %s not found" % id):
            fetch_amplicon(id_, 'deblur', 100)

    def test_study_without_processing_type(self):
        id_ = "10313"  # has 18S and ITS
        proctype = 'deblur'
        trim = 100

        msg = "study %s does not have %s-%dnt 16S data" % (id_, proctype, trim)
        with self.assertRaisesRegex(ValueError, msg):
            fetch_amplicon(id_, proctype, trim)

    def test_study_does_not_have_trim_length(self):
        id_ = "10317"  # has 18S and ITS
        proctype = 'deblur'
        trim = 250

        msg = "study %s does not have %s-%dnt 16S data" % (id_, proctype, trim)
        with self.assertRaisesRegex(ValueError, msg):
            fetch_amplicon(id_, proctype, trim)


if __name__ == '__main__':
    unittest.main()

