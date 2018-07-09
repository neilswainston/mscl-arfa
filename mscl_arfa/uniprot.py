'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import os
from urllib import urlopen, urlretrieve
import xml.sax

import pandas as pd


class UniprotHandler(xml.sax.ContentHandler):
    '''Handler for Uniprot XML files.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.__in_ebml = False
        self.__gen_dna_id = None
        self.__gen_dna_ids = set([])

    def get_gen_dna_ids(self):
        '''Get genomic DNA ids.'''
        return self.__gen_dna_ids

    def startElement(self, name, attrs):
        if name == 'dbReference' and attrs.get('type', None) == 'EMBL':
            self.__in_ebml = True
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'protein sequence ID':
            self.__gen_dna_id = attrs['value']
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'molecule type' \
                and attrs.get('value', None) == 'Genomic_DNA':
            self.__gen_dna_ids.add(self.__gen_dna_id)

    def endElement(self, name):
        if name == 'dbReference':
            self.__in_ebml = False


def get_gen_dna_ids(uniprot_id):
    '''Parse Uniprot XML file.'''
    parser = xml.sax.make_parser()

    handler = UniprotHandler()
    parser.setContentHandler(handler)
    parser.parse(urlopen(
        'https://www.uniprot.org/uniprot/' + uniprot_id + '.xml'))

    return handler.get_gen_dna_ids()


def get_data(name, query, out_dir):
    '''Extract data from Uniprot.'''

    # Download Uniprot data for EC term:
    uniprot_csv = os.path.join(out_dir, name + '_uniprot.tsv')

    if not os.path.exists(uniprot_csv):
        query_str = query + \
            '&format=tab' + \
            '&columns=id,entry name,organism,organism-id'

        url = 'http://www.uniprot.org/uniprot/?query=' + query_str

        urlretrieve(url, uniprot_csv)

    # Read Uniprot data into Dataframe:
    df = pd.read_csv(uniprot_csv, sep='\t')
    df.name = name

    return df
