# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import pandas as pd
import biom
import numpy as np

from skbio import OrdinationResults
from q2_american_gut._visualizer import _insanity_checker


class AGReportTests(unittest.TestCase):
    def setUp(self):
        self.alpha = pd.Series([1, 2, 3], index=list('abc'))

        data = np.asarray([[0, 0, 1], [1, 3, 42]])
        self.biom = biom.Table(data, ['O1', 'O2'], ['a', 'b', 'c'])

        eigvals = [0.51236726, 0.30071909, 0.26791207]
        proportion_explained = [0.2675738328, 0.157044696, 0.1399118638]
        sample_ids = ['a', 'b', 'c']
        axis_labels = ['PC%d' % i for i in range(1, 4)]
        np.random.seed(11)
        data = np.random.randn(3, 3)

        expected_results = OrdinationResults(
            short_method_name='PCoA',
            long_method_name='Principal Coordinate Analysis',
            eigvals=pd.Series(eigvals, index=axis_labels),
            samples=pd.DataFrame(
                data,
                index=sample_ids, columns=axis_labels),
            proportion_explained=pd.Series(proportion_explained,
                                           index=axis_labels))
        self.ordination = expected_results

        self.metadata = pd.DataFrame(data=[[':0', ':)', ':/'],
                                           [':D', 'xD', '<3'],
                                           [';L', ']:->', ':S']],
                                     index=list('abc'),
                                     columns=['foo', 'bar', 'baz'])

    def test_insanity_checks_metadata_error(self):
        pass

    def test_insanity_checks_biom_error(self):
        pass

    def test_insanity_checks_alpha_error(self):
        pass

    def test_insanity_checks_ordination_error(self):
        pass

    def test_sanity_checker(self):
        samples = ['a', 'b', 'c']
        _insanity_checker(samples,
                          self.metadata,
                          self.biom,
                          self.alpha,
                          self.ordination)
        samples = ['a']
        _insanity_checker(samples,
                          self.metadata,
                          self.biom,
                          self.alpha,
                          self.ordination)


if __name__ == "__main__":
    unittest.main()
