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
        
        line = open(str(self)).readline().strip('\n')
        
        return len(line.split('\t')) == 4
        # change this to look through the first five rows and validate that 
        # they are tab delimited and all have the same number of rows


# looking at the source code in model/directory_format, it is unclear to me
# what SingleFileDirectoryFormat does. 
QiitaMetadataDirectoryFormat = model.SingleFileDirectoryFormat(
        'QiitaMetadataDirectoryFormat', 'qiita-metadata.tsv', 
        QiitaMetadataFormat)

#plugin.register_formats(QiitaMetadataFormat, QiitaMetadataDirectoryFormat)
