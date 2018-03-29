# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model

from q2_american_gut.plugin_setup import plugin

class QiitaMetadataFormat(model.TextFileFormat):
    def sniff(self):
        # take the first line of the qiita metadata file
        line = open(str(self)).readline().strip('\n')
        # I need to know if there are some exact specifications that the 
        # metadata file needs to fulfill. In that case we can set up an 
        # expected variable as in q2_quality_filter format Otherwise right 
        # here I am going to check if the first line has 4 parts as in 
        # qiita-metadata.tsv
        return len(line.split('\t')) == 4

# looking at the source code in model/directory_format, it is unclear to me
# what SingleFileDirectoryFormat does. 
QiitaMetadataDirectoryFormat = \
        model.SingleFileDirectoryFormat('QiitaMetadataDirectoryFormat', 
                'qiita-metadata.tsv', QiitaMetadataFormat)
