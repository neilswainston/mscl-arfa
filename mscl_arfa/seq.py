'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
from difflib import SequenceMatcher
import itertools
import os
import sys

from mscl_arfa.ena import get_start_end_comp
from mscl_arfa.uniprot import get_data, get_gen_dna_ids
import numpy as np
import pandas as pd


def run(out_dir):
    '''Run script.'''
    # Make output directory:
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Extract data from Uniprot:
    arfa_query = 'database:(type:pfam id:PF03889)'
    arfa_df = _get_uniprot_data('arfA', arfa_query, out_dir)

    mscl_query = 'database:(type:pfam id:PF01741)'
    mscl_df = _get_uniprot_data('mscL', mscl_query, out_dir)

    df = pd.merge(arfa_df, mscl_df, left_index=True, right_index=True)

    # Get start, end, is complement:
    df = _get_start_ends(df, out_dir)

    # Calculate is overlaps:
    df = _calc_overlaps(df, out_dir)

    df = _pair_genomic_dna_ids(df, out_dir)

    _filter(df, out_dir)


def _get_uniprot_data(name, query, out_dir):
    '''Get Uniprot data.'''

    # Get Uniprot data from API:
    df = get_data(name, query, out_dir)

    # Get genomic DNA id from Uniprot XML:
    gen_dna_csv = os.path.join(out_dir, name + '_gen_dna.csv')

    if not os.path.exists(gen_dna_csv):
        data = []

        for _, row in df.iterrows():
            uniprot_id = row['Entry']

            data.extend([[uniprot_id, prot_seq_id]
                         for prot_seq_id in get_gen_dna_ids(uniprot_id)])

        gen_dna_id_df = pd.DataFrame(data, columns=['Entry', 'genomic_dna_id'])
        gen_dna_id_df.to_csv(gen_dna_csv, encoding='utf8', index=False)
    else:
        gen_dna_id_df = pd.read_csv(gen_dna_csv, encoding='utf8')

    # Merge Uniprot data with genomic id data:
    df = df.merge(gen_dna_id_df, on='Entry')

    # Reset indices and columns:
    df.set_index(['Organism', 'Organism ID'], inplace=True)
    df.columns = pd.MultiIndex.from_tuples([[name, col] for col in df.columns])
    return df


def _get_start_ends(df, out_dir):
    '''Get all start, end, is complement for dataframe.'''
    start_ends_csv = os.path.join(out_dir, 'start_ends.csv')

    if not os.path.exists(start_ends_csv):
        for level in df.columns.levels[0]:
            data = [_get_start_end(genomic_dna_id)
                    for genomic_dna_id in df[level]['genomic_dna_id'].unique()]

            cols = pd.MultiIndex.from_tuples([[level, column]
                                              for column in ['genomic_dna_id',
                                                             'start',
                                                             'end',
                                                             'is_complement']])

            start_end_df = pd.DataFrame(data, columns=cols)

            df = df.reset_index().merge(start_end_df).set_index(df.index.names)

        df.to_csv(start_ends_csv, encoding='utf-8')
    else:
        df = pd.read_csv(start_ends_csv, index_col=[0, 1], header=[0, 1])

    return df


def _get_start_end(genomic_dna_id):
    '''Get start, end, is complement from genomic DNA id.'''
    return [genomic_dna_id] + list(get_start_end_comp(genomic_dna_id))


def _calc_overlaps(df, out_dir):
    '''Calculate overlaps.'''
    overlaps_csv = os.path.join(out_dir, 'overlaps.csv')

    if not os.path.exists(overlaps_csv):
        is_overlaps = []

        levels = dict(zip(df.columns.get_level_values(0),
                          df.columns.get_level_values(1))).keys()

        for _, row in df.iterrows():
            vals = [[row[level]['start'],
                     row[level]['end'],
                     row[level]['is_complement']]
                    for level in levels]
            is_overlaps.append(_calc_overlap(vals[0], vals[1]))

        df['common', 'overlap'] = is_overlaps
        df.to_csv(overlaps_csv, encoding='utf-8')
    else:
        df = pd.read_csv(overlaps_csv, index_col=[0, 1], header=[0, 1])

    return df


def _calc_overlap(left, right):
    '''Calculate overlap if genes are on complementary strands.'''
    if all(pd.notna(left)) and all(pd.notna(right)) and left[2] ^ right[2]:
        range_left = range(int(left[0]), int(left[1]))
        range_right = range(int(right[0]), int(right[1]))
        intersection = set(range_left).intersection(range_right)
        return len(intersection)

    return 0


def _pair_genomic_dna_ids(df, out_dir):
    paired_gen_ids_csv = os.path.join(out_dir, 'paired_gen_ids.csv')

    if not os.path.exists(paired_gen_ids_csv):
        '''Pair genomic_dna_ids.'''
        df['common', 'gen_data_id_sim'] = \
            df.apply(__score_gen_data_id_similarity, axis=1)
        df.to_csv(paired_gen_ids_csv, encoding='utf-8')
    else:
        df = pd.read_csv(paired_gen_ids_csv, index_col=[0, 1], header=[0, 1])

    return df


def __score_gen_data_id_similarity(row):
    '''Score genomic data id similarity.'''
    gen_data_ids = list(row.xs('genomic_dna_id', level=1))

    similarities = [SequenceMatcher(None, id_1, id_2).ratio()
                    for id_1, id_2 in itertools.combinations(gen_data_ids, 2)]

    return np.mean(similarities)


def _filter(df, out_dir):
    '''Filter data based on overlap > 0 and paired genomic DNA ids,
    retaining "best" pairs per organism.'''
    filtered_csv = os.path.join(out_dir, 'filtered.csv')

    if not os.path.exists(filtered_csv):
        series = []

        filtered_df = df[df['common', 'overlap'] > 0]

        for _, pairs_df in filtered_df.groupby(filtered_df.index):
            series_df = \
                pairs_df.sort_values(by=[('common', 'gen_data_id_sim')],
                                     ascending=False)

            series.append(series_df.iloc[0])

        filtered_df = \
            pd.DataFrame(series).sort_values(by=[('common',
                                                  'gen_data_id_sim')],
                                             ascending=False)

        filtered_df.to_csv(filtered_csv, encoding='utf-8')
    else:
        filtered_df = pd.read_csv(filtered_csv,
                                  index_col=[0, 1], header=[0, 1])

    return filtered_df


def main(args):
    '''main method.'''
    run(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
