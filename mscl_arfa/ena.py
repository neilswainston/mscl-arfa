'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import re
import sys
import xml.sax

_RE = r'(\w+)\(.*:(\d+)..(\d+)\)'


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
            groups = re.match(_RE, location)
            self.__complement = groups.group(1) == 'complement'
            self.__start = int(groups.group(2))
            self.__end = int(groups.group(3))


def parse(filename):
    '''Parse ENA XML file.'''
    parser = xml.sax.make_parser()

    handler = EnaHandler()
    parser.setContentHandler(handler)

    with open(filename) as fle:
        parser.parse(fle)

    return handler.get_start(), handler.get_end(), handler.is_complement()


def main(args):
    '''main method.'''
    print parse(args[0])


if __name__ == '__main__':
    main(sys.argv[1:])
