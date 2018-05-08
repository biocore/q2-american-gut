# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2
from q2_american_gut.plugin_setup import plugin
from q2_american_gut import QiitaMetadataFormat, QiitaMetadata


@plugin.register_transformer
def _1(data: pd.DataFrame) -> QiitaMetadataFormat:
    ff = QiitaMetadataFormat()
    with ff.open() as fh:
        data.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _2(ff: QiitaMetadataFormat) -> pd.DataFrame:
    with ff.open() as fh:
        df = pd.read_csv(fh, sep='\t', dtype=object)
        return df


@plugin.register_transformer
def _3(ff: QiitaMetadataFormat) -> qiime2.Metadata:
    with ff.open() as fh:
        df = pd.read_csv(fh, sep='\t', dtype='object')
        df.set_index('#SampleID', inplace=True)
        return qiime2.Metadata(df)


@plugin.register_transformer
def _4(data: qiime2.Metadata) -> QiitaMetadataFormat:
    ff = QiitaMetadataFormat()
    md_df = data.to_dataframe()
    with ff.open() as fh:
        md_df.to_csv(fh, sep='\t', header=True)
    return ff


@plugin.register_transformer
def _5(data: qiime2.Metadata) -> QiitaMetadata:
    pass


@plugin.register_transformer
def _6(data: QiitaMetadata) -> qiime2.Metadata:
    pass
