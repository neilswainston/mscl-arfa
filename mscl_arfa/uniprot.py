'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=protected-access
import os
import ssl
from urllib.request import urlopen, urlretrieve
import xml.sax

import pandas as pd


# Apparently a dirty security hack:
ssl._create_default_https_context = ssl._create_unverified_context


class UniprotHandler(xml.sax.ContentHandler):
    '''Handler for Uniprot XML files.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.__in_ebml = False
        self.__gen_dna_id = None
        self.__embl_id = None
        self.__gen_dna_ids = set([])
        self.__embl_ids = {}

    def get_gen_dna_ids(self):
        '''Get genomic DNA ids.'''
        return self.__gen_dna_ids

    def get_embl_ids(self):
        '''Get EMBL ids.'''
        return self.__embl_ids

    def startElement(self, name, attrs):
        if name == 'dbReference' and attrs.get('type', None) == 'EMBL':
            self.__in_ebml = True
            self.__embl_id = attrs['id']
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'protein sequence ID':
            self.__gen_dna_id = attrs['value']
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'molecule type' \
                and attrs.get('value', None) == 'Genomic_DNA':
            self.__gen_dna_ids.add(self.__gen_dna_id)
            self.__embl_ids[self.__embl_id] = self.__gen_dna_id

    def endElement(self, name):
        if name == 'dbReference':
            self.__in_ebml = False


def get_gen_dna_ids(uniprot_id):
    '''Get DNA ids.'''
    return _parse(uniprot_id).get_gen_dna_ids()


def get_embl_ids(uniprot_id):
    '''Get EMBL ids.'''
    return _parse(uniprot_id).get_embl_ids()


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


def _parse(uniprot_id):
    '''Parse xml.'''
    parser = xml.sax.make_parser()
    handler = UniprotHandler()
    parser.setContentHandler(handler)
    fle = urlopen('https://www.uniprot.org/uniprot/' + uniprot_id + '.xml')
    parser.parse(fle)
    return handler
