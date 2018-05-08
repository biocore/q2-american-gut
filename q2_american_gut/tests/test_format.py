
# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import shutil

from q2_american_gut.plugin_setup import QiitaMetadataFormat, \
                                         QiitaMetadataDirectoryFormat
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError

class TestFormats(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_qiita_metadata_validate_positive(self):
        filepath = self.get_data_path('qiita-metadata.tsv')
        format = QiitaMetadataFormat(filepath, mode='r')

        format.validate()

    def test_qiita_metadata_format_validate_negative(self):
        filepath = self.get_data_path('not-qiita-metadata.tsv')
        format = QiitaMetadataFormat(filepath, mode='r')

        with self.assertRaisesRegex(ValidationError, 'QiitaMetadataFormat'):
            format.validate()
    
    def test_qiita_metadata_dir_fmt_validate_positive(self):
        filepath = self.get_data_path('qiita-metadata.tsv')
        shutil.copy(filepath, self.temp_dir.name)
        format = QiitaMetadataDirectoryFormat(self.temp_dir.name, mode='r')

        format.validate()


if __name__ == '__main__':
    unittest.main()
