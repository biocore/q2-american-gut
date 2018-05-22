# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_american_gut._type import QiitaMetadata
from q2_american_gut._format import QiitaMetadataDirectoryFormat
from qiime2.plugin.testing import TestPluginBase


class TestTypes(TestPluginBase):
    package = "q2_american_gut.tests"

    def test_qiita_metadata_semantic_type_registration(self):
        self.assertRegisteredSemanticType(QiitaMetadata)

    def test_qiita_metadata_semantic_type_directory_format(self):
        self.assertSemanticTypeRegisteredToFormat(
            QiitaMetadata, QiitaMetadataDirectoryFormat)


if __name__ == '__main__':
    unittest.main()
