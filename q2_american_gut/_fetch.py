# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import skbio


def fetch_amplicon(qiita_study_id: int, processing_type: str, trim_length: int,
                   threads: int=1, debug: bool=False) -> (biom.Table,
                                                          pd.DataFrame,
#                                                          pd.DataFrame,
                                                          skbio.TreeNode):
    pass
