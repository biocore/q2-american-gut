# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import SemanticType
#from q2_american_gut.plugin_setup import plugin
# I am not sure here if this type needs other parameters in the SemTyp 
# function. Quality filters uses only this but sample data is a little more
# complicated
QiitaMetadata = SemanticType('QiitaMetadata')

#plugin.register_semantic_types(QiitaMetadata)


