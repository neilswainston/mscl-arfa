'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import re
import sys
from urllib import urlopen
import xml.sax
_RE = r'(\d+)\.\.(\d+)'


class EnaHandler(xml.sax.ContentHandler):
    '''Handler for ENA XML files.'''

    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.__start = None
        self.__end = None
        self.__complement = None

    def get_start(self):
        '''Get start.'''
        return self.__start

    def get_end(self):
        '''Get end.'''
        return self.__end

    def is_complement(self):
        '''Get complement flah.'''
        return self.__complement

    def startElement(self, name, attrs):
        if name == 'feature' and attrs.get('name', None) == 'CDS':
            location = attrs['location']
            self.__complement = location.startswith('complement')
            matches = re.findall(_RE, location)
            self.__start = int(matches[0][0])
            self.__end = int(matches[0][1])


def get_start_end_comp(ena_id):
    '''Get start, end, is complement.'''
    url = 'https://www.ebi.ac.uk/ena/data/view/%s&display=xml&download=xml' \
        % ena_id
    return _parse(urlopen(url))


def parse(filename):
    '''Parse ENA XML file.'''
    with open(filename) as fle:
        return _parse(fle)


def _parse(src):
    '''Parse source data.'''
    parser = xml.sax.make_parser()
    handler = EnaHandler()
    parser.setContentHandler(handler)
    parser.parse(src)
    return handler.get_start(), handler.get_end(), handler.is_complement()


def main(args):
    '''main method.'''
    print parse(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
