'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=broad-except
# pylint: disable=wrong-import-order
import sys
import tempfile

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import requests

from mscl_arfa import uniprot
import pandas as pd


def get_seqs(in_filename, extension, out_filename):
    '''Get mscl sequences.'''
    in_df = pd.read_csv(in_filename, header=[0, 1])
    uniprot_ids = in_df[('mscL', 'Entry name')].dropna().unique()

    records = []

    for uniprot_id in uniprot_ids:
        print('Extracting data for %s' % uniprot_id)

        try:
            for entry_id, prot_id in uniprot.get_embl_ids(uniprot_id).items():
                seq = _get_genbank(entry_id, prot_id, extension)

                if seq:
                    records.append(SeqRecord(seq, id=entry_id, name=entry_id,
                                             description=entry_id))
        except Exception:
            print('Unable to extract sequence data for %s' % uniprot_id)

    with open(out_filename, 'w') as fle:
        SeqIO.write(records, fle, 'fasta')


def _get_genbank(entry_id, prot_id, extension, db_id='nuccore'):
    '''Get Genbank data.'''
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + \
        'db=%s&id=%s&rettype=gb&retmode=text' % (db_id, entry_id)

    tmpfile = tempfile.NamedTemporaryFile(delete='False')
    resp = requests.get(url)

    with open(tmpfile.name, 'w') as fle:
        fle.write(resp.text)

    genbank = SeqIO.read(tmpfile.name, 'genbank')
    seq = genbank.seq

    for feature in [feat for feat in genbank.features
                    if feat.type == 'CDS'
                    and prot_id in feat.qualifiers.get('protein_id', [])]:
        loc = feature.location
        if loc.strand == 1:
            return seq[loc.start:loc.end + extension]
        # else:
        return seq[loc.start - extension:loc.end].reverse_complement()

    return None


def main(args):
    '''main method.'''
    get_seqs(args[0], int(args[1]), args[2])


if __name__ == '__main__':
    main(sys.argv[1:])
