# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pkg_resources import Requirement, resource_filename

import qiime2
from q2_types.feature_data import DNAIterator
import biom
import pandas as pd
import skbio
from redbiom import search, fetch, summarize
import redbiom.search
import redbiom.summarize


CLASSIFIER = (Requirement.parse('q2_american_gut'),
              'q2_american_gut/assets/gg-13-8-99-515-806-nb-classifier.qza')


def _determine_context(processing_type, trim_length):
    """Determine the processing context in redbiom to use

    Parameters
    ----------
    processing_type : str, {deblur, closed-reference}
        The sequence processing type to look for.
    trim_length : int
        The fragment length to look for.

    Notes
    -----
    This method is restricted to 16S V4 data.

    Raises
    ------
    ValueError
        If a suitable context could not be found.
    """
    contexts = redbiom.summarize.contexts()

    found = None
    for context in contexts.ContextName:
        context = context.lower()

        # american gut amplicon data are only 16S v4
        if '16s' not in context or 'v4' not in context:
            continue

        # TODO: the release candidate version of the qiita redbiom database
        # has updated context names which include the reference used.
        # if we're closed reference, restrict to greengenes as a reference
        #if processing_type != 'deblur' and 'greengenes' not in context:
        #    continue

        if processing_type in context and "%dnt" % trim_length in context:
            found = context

    if found is None:
        msg = "Cannot find %s-%dnt for 16S data" % (processing_type, trim_length)
        raise ValueError(msg)

    return found


def _get_featuredata_from_table(table):
    """Extract the observations and interpret as skbio.DNA"""
    if table.is_empty():
        raise ValueError("No features")

    return DNAIterator((skbio.DNA(i, metadata={'id': i})
                                  for i in table.ids(axis='observation')))


def _assign_taxonomy(table, threads):
    from qiime2.plugins import feature_classifier
    dna_iter = _get_featuredata_from_table(table)

    classifier = Artifact.load(resource_filename(*CLASSIFIER))
    tax, = feature_classifier.methods.classify_sklearn(dna_iter, classifier,
                                                       n_jobs=threads)
    return tax

def _get_taxonomy(table):
    #might need test, test for a proper biom table
    pass

def _get_closed_reference_phylogeny():

    pass

def _insert_fragments(table, threads):
    from qiime2.plugins import fragment_insertion
    dna_iter = _get_featuredata_from_table(table)
    return skbio.TreeNode()


def fetch_amplicon(qiita_study_id: str, processing_type: str, trim_length: int,
                   threads: int=1, debug: bool=False) -> (biom.Table,
                                                          pd.DataFrame,
                                                          pd.DataFrame,
                                                          skbio.TreeNode):

    if not qiita_study_id.isdigit():
        raise ValueError("ID %s is not a qiita ID" % qiita_study_id)

    if int(qiita_study_id) < 1:
        raise ValueError("ID %s is not a qiita ID" % qiita_study_id)

    context = _determine_context(processing_type, trim_length)

    query = "where qiita_study_id==%s" % qiita_study_id
    samples = redbiom.search.metadata_full(query)

    if not samples:
        raise ValueError("study ID %s has no samples or is not a Qiita ID")

    table, ambiguity_map_tab = redbiom.fetch.data_from_samples(context, samples)
    md, ambiguity_map_md = redbiom.fetch.metadata_from_samples(context, samples)

    # TODO: do not execute these steps if we're pulling out closed reference
    # as the taxonomy and tree already exist
    if processing_type == 'closed-reference':
        taxonomy = _get_taxonomy(table)
        phylogeny = _get_closed_reference_phylogeny()
        pass
    else:
        taxonomy = _assign_taxonomy(table, threads) # for closed ref you do not
        phylogeny = _insert_fragments(table, threads) # assign tax or insert frag

    # qiita study 2136

    return table, taxonomy, md, phylogeny

