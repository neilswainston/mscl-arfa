'''
SYNBIOCHEM (c) University of Manchester 2018

SYNBIOCHEM is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import unittest

from mscl_arfa.uniprot import get_prot_seq_ids


class Test(unittest.TestCase):
    '''Test class for parser.'''

    def test_get_prot_seq_ids(self):
        '''Tests get_prot_seq_ids method.'''
        expected = [u'AAC50125.1', u'BAA23527.1', u'AAP88775.1', u'EAW65858.1',
                    u'EAW65856.1', u'EAW65859.1']
        self.assertEqual(get_prot_seq_ids('P43699'), expected)
