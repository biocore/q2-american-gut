# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd

from q2_american_gut._visualizer import _insanity_checker


class AGReportTests(unittest.TestCase):
    def setUp(self):
        self.alpha = None
        self.biom = None
        self.ordination = None
        self.metadata = None

    def test_insanity_checks_metadata_error(self):
        pass

    def test_insanity_checks_biom_error(self):
        pass

    def test_insanity_checks_alpha_error(self):
        pass

    def test_insanity_checks_ordination_error(self):
        pass

if __name__ == "__main__":
    unittest.main()
