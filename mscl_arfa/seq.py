'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import os
import sys

from mscl_arfa.ena import get_start_end_comp
from mscl_arfa.uniprot import get_data, get_gen_dna_ids
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
    _is_overlaps(df)

    df.to_csv(os.path.join(out_dir, 'out.csv'), encoding='utf-8')


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

        df.to_csv(start_ends_csv)
    else:
        df = pd.read_csv(start_ends_csv, index_col=[0, 1], header=[0, 1])

    return df


def _get_start_end(genomic_dna_id):
    '''Get start, end, is complement from genomic DNA id.'''
    return [genomic_dna_id] + list(get_start_end_comp(genomic_dna_id))


def _is_overlaps(df):
    '''Get is overlap data.'''
    is_overlaps = []

    for _, row in df.iterrows():
        vals = [[row[level]['start'],
                 row[level]['end'],
                 row[level]['is_complement']]
                for level in dict(zip(df.columns.get_level_values(0),
                                      df.columns.get_level_values(1))).keys()]
        is_overlaps.append(_is_overlap(vals[0], vals[1]))

    df['overlap'] = is_overlaps


def _is_overlap(left, right):
    '''Calculate whether genes are complementary and overlap.'''
    if all(pd.notna(left)) and all(pd.notna(right)) and left[2] ^ right[2]:
        range_left = range(int(left[0]), int(left[1]))
        range_right = range(int(right[0]), int(right[1]))
        intersection = set(range_left).intersection(range_right)
        return len(intersection)

    return 0


def main(args):
    '''main method.'''
    run(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
