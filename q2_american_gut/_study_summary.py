# ----------------------------------------------------------------------------
# Copyright (c) 2012-2018, American Gut Project development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_american_gut._fetch import fetch_amplicon
import skbio
from bloom import remove_seqs, trim_seqs


def summarize_study(ctx, qiita_study_id, processing_type, trim_length,
                   threads=1, debug=False):

    """Summarizes data from the input qiita study id

    Parameters
    ----------
    ctx: qiime2.sdk.Context
        The QIIME2 (note: NOT redbiom) context; is required first parameter of
        any function that defines a QIIME2 pipeline
    qiita_study_id: str
        A string containing a qiita study id; should be a digit >= 1
    processing_type : {'deblur', 'closed-reference'}
        The processing type which determines how we resolve a phylogeny.
    trim_length : int
        The fragment length to look for.
    threads : int
        The number of threads to use
    debug : bool, optional
        If debug, use a small tree for insertion to avoid spinup
        time.

    Returns
    -------
    TODO

    Notes
    -----
    This method is restricted to 16S V4 data.
    This method selects only data from the illumina platform.
    """

    # minimum acceptable target: being able to provide PCoA of unweighted
    # unifrac distances that someone can interact with to explore study
    # metadata

    # For input qiita study id:
    # 1. Remove blooms
    # 2. Rarefy table to some reasonably permissive depth
    # 3. Compute unweighted unifrac over full feature table
    # 4. Produce PCoA
    # Additional:
    # 5. Specify initial category paint of PCoA plot
    # 6. Produce basic statistics
    #   like # of samples, sequences per sample, break down by a parameterized
    #   metadata category (e.g., body site)
    # 7. Time stamp the production of results

    # Integrating participants data with the human microbiome project:
    # need to do that through closed-reference.
    # Don't remove blooms and hope it doesn't make a difference.
    # If not, filter out exact blooms, then do closed-reference picking to make
    # comparable with HMP.

    # table is a FeatureTable[Frequency]
    # taxonomy is a FeatureData[Taxonomy]
    # mtadata is a QiitaMetadata
    # phylogeny is a Phylogeny[Rooted]
    table, taxonomy, mtadta, phylogeny = fetch_amplicon(
        ctx, qiita_study_id, processing_type, trim_length, threads=threads,
        debug=debug)

    # TODO: How to automatically pick a "reasonably permissive sampling depth"?
    sampling_depth = None  # NB: not a real value, won't work

    # mostly based on https://github.com/qiime2/q2-diversity/q2_diversity/
    # _core_metrics.py#L45-L47
    rarefy = ctx.get_action('feature_table', 'rarefy')
    beta_phylogenetic = ctx.get_action('diversity', 'beta_phylogenetic')
    pcoa = ctx.get_action('diversity', 'pcoa')
    emperor_plot = ctx.get_action('emperor', 'plot')

    rarefied_table, = rarefy(table=table, sampling_depth=sampling_depth)

    # below based heavily on https://github.com/qiime2/q2-diversity/
    # q2_diversity/_core_metrics.py#L57-L69
    distance_matrix_list = beta_phylogenetic(
        table=rarefied_table, phylogeny=phylogeny, metric='unweighted_unifrac',
        n_jobs=threads)
    unw_unifrac_dm = distance_matrix_list[0]

    pcoas = pcoa(distance_matrix=unw_unifrac_dm)
    unw_unifrac_pcoa = pcoas[0]

    # TODO: specify custom_axes=<a List[Str]> param
    unw_unifrac_emp_plot = emperor_plot(pcoa=unw_unifrac_pcoa, metadata=mtadta)

    return (
        table, taxonomy, mtadta, phylogeny, rarefied_table, unw_unifrac_dm,
        unw_unifrac_pcoa, unw_unifrac_emp_plot)


def _remove_blooms(biom_table, trim_length):
    # 1. Remove blooms
    #   See knightlab-analyses/bloom-analyses repo: see data/newbloom.all.fna:
    #   these are 150 nt sequences identified as sequences associated with
    #   organisms that seem to bloom during shipping.
    #   Literally delete these rows of the biom table and move forward.
    #   NB: Since these are 150 nt, if filtering at 90 nt data, need to trim
    #   bloom sequences to 90 nt so can determine if have exact match
    #   (yes, doing exact match).

    # TODO: get the location of the file of bloom fasta sequences,
    # newbloom.all.fna
    bloom_seqs_fp = None  # Not a real value, won't work

    # below based on https://github.com/knightlab-analyses/bloom-analyses/blob/
    # master/ipynb/bloom_example.ipynb
    # Theoretically, I could do this by calling
    # bloom.filter_seqs_from_biom.main(); however (a) a little bit of
    # refactoring of that code would be necessary to allow the main method
    # to take parameters via a code call rather than from sys.argv
    # (see https://stackoverflow.com/a/14500228) and (b) that method expects
    # files as its biom inputs and outputs, which may be overkill here as
    # we already have the stuff in memory.

    bloom_seqs = skbio.read(bloom_seqs_fp, format='fasta')
    trimmed_bloom_seqs = trim_seqs(bloom_seqs, seqlength=trim_length)
    bloom_removed_biom_table = remove_seqs(biom_table, trimmed_bloom_seqs)
    return bloom_removed_biom_table
