'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from urllib import urlopen
import xml.sax


class UniprotHandler(xml.sax.ContentHandler):
    '''Handler for Uniprot XML files.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.__in_ebml = False
        self.__prot_seq_id = None
        self.__prot_seq_ids = []

    def get_prot_seq_ids(self):
        '''Get protein sequence ids.'''
        return self.__prot_seq_ids

    def startElement(self, name, attrs):
        if name == 'dbReference' and attrs.get('type', None) == 'EMBL':
            self.__in_ebml = True
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'protein sequence ID':
            self.__value = attrs['value']
        elif self.__in_ebml and name == 'property' \
                and attrs.get('type', None) == 'molecule type' \
                and attrs.get('value', None) == 'Genomic_DNA':
            self.__prot_seq_ids.append(self.__value)

    def endElement(self, name):
        if name == 'dbReference':
            self.__in_ebml = False


def get_prot_seq_ids(uniprot_id):
    '''Parse ENA XML file.'''
    parser = xml.sax.make_parser()

    handler = UniprotHandler()
    parser.setContentHandler(handler)
    parser.parse(urlopen(
        'https://www.uniprot.org/uniprot/' + uniprot_id + '.xml'))

    return handler.get_prot_seq_ids()


def main(args):
    '''main method.'''
    print get_prot_seq_ids(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
