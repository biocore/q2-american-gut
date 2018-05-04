# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from pkg_resources import Requirement, resource_filename

import qiime2
from q2_types.feature_data import (DNAIterator,
                                   AlignedDNASequencesDirectoryFormat)
from q2_types.tree import NewickFormat
import biom
import pandas as pd
import skbio
import redbiom.search
import redbiom.summarize
import redbiom.fetch

from q2_american_gut._type import QiitaMetadata

CLASSIFIER = (Requirement.parse('q2_american_gut'),
              'q2_american_gut/assets/gg-13-8-99-515-806-nb-classifier.qza')
GG_TREE = (Requirement.parse('q2_american_gut'),
          'q2_american_gut/assets/97_otus.tree')
DEBUG_PHY = (Requirement.parse('q2_american_gut'),
             'q2_american_gut/assets/reference_phylogeny_tiny.qza')
DEBUG_ALN = (Requirement.parse('q2_american_gut'),
            'q2_american_gut/assets/reference_alignment_tiny.qza')


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
        test_ctx = context.lower()

        # american gut amplicon data are only 16S v4
        if '16s' not in test_ctx or 'v4' not in test_ctx:
            continue

        # TODO: the release candidate version of the qiita redbiom database
        # has updated context names which include the reference used.
        # if we're closed reference, restrict to greengenes as a reference
        #if processing_type != 'deblur' and 'greengenes' not in test_ctx:
        #    continue

        if processing_type in test_ctx and "%dnt" % trim_length in test_ctx:
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


def _fetch_taxonomy(processing_type, table, threads):
    """Fetch taxonomy based on the processing type

    Parameters
    ----------
    processing_type : {'deblur', 'closed-reference'}
        The processing type which determines how we obtain taxonimic detail.
    table : biom.Table
        The FrequencyTable
    threads : int
        The number of threads to use

    Returns
    -------
    FeatureData[Taxonomy]
    """
    if processing_type == 'deblur':
        
        from qiime2.plugins import feature_classifier
        dna_iter = _get_featuredata_from_table(table)

        reads = qiime2.Artifact.import_data('FeatureData[Sequence]', dna_iter)
        classifier = qiime2.Artifact.load(resource_filename(*CLASSIFIER))
        tax, = feature_classifier.methods.classify_sklearn(reads,
                                                           classifier,
                                                           n_jobs=threads)
    else:
        tax_table = pd.DataFrame([(i, m['taxonomy'], 1.0)
                                  for v, i, m in table.iter(axis='observation',
                                                            dense=False)],
                                 columns=['Feature ID', 'Taxon', 'Confidence'])
        tax_table.set_index('Feature ID', inplace=True)
        tax = qiime2.Artifact.import_data('FeatureData[Taxonomy]', tax_table)

    return tax

def _fetch_phylogeny(processing_type, table, threads, debug):
    """Fetch phylogeny based on the processing type

    Parameters
    ----------
    processing_type : {'deblur', 'closed-reference'}
        The processing type which determines how we resolve a phylogeny.
    table : biom.Table
        The FrequencyTable
    threads : int
        The number of threads to use
    debug : bool, optional
        If debug, use a small small tree for insertion to avoid spinup
        time.

    Returns
    -------
    Phylogeny[Rooted]
    """
    if processing_type == 'deblur':
        from qiime2.plugins.fragment_insertion.methods import sepp

        # insertion is extremely expensive to spin up (20min) so short
        # circuit if debug
        if debug:
            ref_phy = qiime2.Artifact.load(resource_filename(*DEBUG_PHY))
            #ref_phy = ref_phy.view(NewickFormat)

            ref_aln = qiime2.Artifact.load(resource_filename(*DEBUG_ALN))
            #ref_aln = ref_aln.view(AlignedDNASequencesDirectoryFormat)
        else:
            ref_phy = None
            ref_aln = None

        dna_iter = _get_featuredata_from_table(table)
        reads = qiime2.Artifact.import_data('FeatureData[Sequence]', dna_iter)
        tree, _ = sepp(reads, threads=threads,
                       reference_alignment=ref_aln,
                       reference_phylogeny=ref_phy)

    else:
        tree = skbio.TreeNode.read(resource_filename(*GG_TREE))
        tree = qiime2.Artifact.import_data('Phylogeny[Rooted]', tree)

    return tree


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

    if debug:
        # just take top 10 but be consistent about which 10 are obtained
        samples = set(sorted(samples)[:10])

    table, ambiguity_map_tab = redbiom.fetch.data_from_samples(context, samples)
    md, ambiguity_map_md = redbiom.fetch.sample_metadata(samples, context=context)
    md.set_index('#SampleID', inplace=True)

    taxonomy = _fetch_taxonomy(processing_type, table, threads)
    phylogeny = _fetch_phylogeny(processing_type, table, threads, debug)

    md = qiime2.Artifact.import_data(QiitaMetadata, md)
    md = qiime2.Metadata.from_artifact(md)
    table = qiime2.Artifact.import_data('FeatureTable[Frequency]', table)

    return table, taxonomy, md, phylogeny

