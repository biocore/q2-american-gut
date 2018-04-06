# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model

class QiitaMetadataFormat(model.TextFileFormat):
    def sniff(self):
        lines = []
        fn = open(str(self))
        for i in range(5):
            lines.append(fn.readline())
        length = len(lines[0].split('\t'))
        if length == 1:
            return False
        for line in lines:
            if len(line.split('\t')) != length:
                return False
        return True

QiitaMetadataDirectoryFormat = model.SingleFileDirectoryFormat(
        'QiitaMetadataDirectoryFormat', 'qiita-metadata.tsv', 
        QiitaMetadataFormat)
