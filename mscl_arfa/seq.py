'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
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
    arfa_df = _get_data('arfA', arfa_query, out_dir)

    mscl_query = 'database:(type:pfam id:PF01741)'
    mscl_df = _get_data('mscL', mscl_query, out_dir)

    df = pd.merge(arfa_df, mscl_df, left_index=True, right_index=True)
    df.dropna(inplace=True)
    df.to_csv(os.path.join(out_dir, 'out.csv'), encoding='utf-8')


def _get_data(name, query, out_dir):
    '''Get data.'''

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

    # Get start, end, is complement:
    df = df.merge(df['genomic_dna_id'].apply(_get_start_end),
                  left_index=True, right_index=True)

    # Reset indices and columns:
    df.set_index(['Organism', 'Organism ID'], inplace=True)
    df.columns = pd.MultiIndex.from_tuples([[name, col] for col in df.columns])
    return df


def _get_start_end(genomic_dna_id):
    '''Get start, end, is complement from genomic DNA id.'''
    return pd.Series(get_start_end_comp(genomic_dna_id),
                     index=['start', 'end', 'is_complement'])


def _is_overlap(left, right):
    '''Calculate whether genes overlap.'''
    return left[2] ^ right[2] and \
        set(range(left[0], left[1])).intersection(range(right[0], right[1]))


def main(args):
    '''main method.'''
    run(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
