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
import redboiom
from redbiom import search, fetch, summarize


def determine_context(processing_type, trim_len):
    
    pass

def fetch_amplicon(qiita_study_id: str, processing_type: str, trim_length: int,
                   threads: int=1, debug: bool=False) -> (biom.Table,
                                                          pd.DataFrame,
                                                          pd.DataFrame,
                                                          skbio.TreeNode):

    
    if not qiita_study_id.is_digit():
        raise ValueError("Qiita study id is not valid")

    

    
    pass

'''
mission
find a public study in qiita, like one of the gold star ones
pull out the deblur 90 nucleotide data for it
'''
