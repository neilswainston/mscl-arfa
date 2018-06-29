'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import shutil
import sys
import urllib
import pandas as pd


def run(out_dir):
    '''Run script.'''
    # Make output directory:
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    os.makedirs(out_dir)

    # Extract data from Uniprot:
    arfa_query = 'database:(type:pfam id:PF03889)'
    arfa_df = _uniprot_extract('arfA', arfa_query, out_dir)

    mscl_query = 'database:(type:pfam id:PF01741)'
    mscl_df = _uniprot_extract('mscL', mscl_query, out_dir)

    df = pd.merge(arfa_df, mscl_df, left_index=True, right_index=True)
    df.dropna(inplace=True)
    df.to_csv(os.path.join(out_dir, 'out.csv'), encoding='utf-8')


def _uniprot_extract(name, query, out_dir):
    '''Extract data from Uniprot.'''

    # Download Uniprot data for EC term:
    uniprot_csv = os.path.join(out_dir, name + '_uniprot.tsv')

    query_str = query + \
        '&format=tab' + \
        '&columns=entry name,organism,organism-id,' + \
        'database(EMBL)'

    url = 'http://www.uniprot.org/uniprot/?query=' + query_str

    urllib.urlretrieve(url, uniprot_csv)

    # Read Uniprot data into Dataframe:
    df = pd.read_csv(uniprot_csv, sep='\t')
    df.name = name
    df.set_index(['Organism', 'Organism ID'], inplace=True)
    df.columns = pd.MultiIndex.from_tuples([[name, col] for col in df.columns])

    return df


def main(args):
    '''main method.'''
    run(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
