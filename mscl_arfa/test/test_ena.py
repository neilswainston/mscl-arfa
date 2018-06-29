'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import unittest

from mscl_arfa.ena import parse


class Test(unittest.TestCase):
    '''Test class for parser.'''

    def test_parse(self):
        '''Tests parse method.'''
        directory = os.path.dirname(os.path.realpath(__file__))
        xml_filepath = os.path.join(directory, '../../data/AAC21905.xml.txt')

        start, end, complement = parse(xml_filepath)
        self.assertEqual(start, 265624)
        self.assertEqual(end, 265833)
        self.assertTrue(complement)
