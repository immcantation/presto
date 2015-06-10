"""
Unit tests for Annotation module
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'

# Imports
import time
import unittest
from collections import OrderedDict

# Presto imports
from presto.Annotation import collapseAnnotation, mergeAnnotation, renameAnnotation


class TestAnnotation(unittest.TestCase):
    def setUp(self):
        print '-> %s()' % self._testMethodName
        # Annotation dictionaries
        self.ann_dict_1 = OrderedDict([('ID', 'SEQ1'), ('TEST1', 'A,B'), ('TEST2', [1, 2])])
        self.ann_dict_2 = OrderedDict([('ID', 'SEQ2'), ('TEST1', 'C,C'), ('TEST2', 3)])

        # Start clock
        self.start = time.time()

    def tearDown(self):
        # End clock
        t = time.time() - self.start
        print '<- %s() %.3f' % (self._testMethodName, t)

    #@unittest.skip('-> mergeAnnotation() skipped\n')
    def test_mergeAnnotation(self):
        result = mergeAnnotation(self.ann_dict_1, self.ann_dict_2)
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A,B,C,C', result['TEST1'])
        self.assertEqual('1,2,3', result['TEST2'])

        result = mergeAnnotation(self.ann_dict_1, self.ann_dict_2, prepend=True)
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('C,C,A,B', result['TEST1'])
        self.assertEqual('3,1,2', result['TEST2'])

    #@unittest.skip('-> renameAnnotation() skipped\n')
    def test_renameAnnotation(self):
        result = renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED1')
        print result
        self.assertEqual('A,B', result['RENAMED1'])
        self.assertEqual(None, result.get('TEST1', None))

        result = renameAnnotation(self.ann_dict_1, 'TEST1', 'RENAMED2')
        print result
        self.assertEqual('A,B', result['RENAMED2'])
        self.assertEqual(None, result.get('TEST1', None))

    #@unittest.skip('-> collapseAnnotation() skipped\n')
    def test_collapseAnnotation(self):
        result = collapseAnnotation(self.ann_dict_1, 'first')
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('A', result['TEST1'])
        self.assertEqual(1, result['TEST2'])

        result = collapseAnnotation(self.ann_dict_1, 'last')
        print result
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('B', result['TEST1'])
        self.assertEqual(2, result['TEST2'])

        result = collapseAnnotation(self.ann_dict_1, 'set')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual(['A' ,'B'], result['TEST1'])
        self.assertEqual([1, 2], result['TEST2'])
        print result

        result = collapseAnnotation(self.ann_dict_1, 'cat')
        self.assertEqual('SEQ1', result['ID'])
        self.assertEqual('AB', result['TEST1'])
        self.assertEqual('12', result['TEST2'])
        print result

        # TODO:  may be better for min/max/sum to return float/int
        result = collapseAnnotation(self.ann_dict_1, 'min', fields='TEST2')
        self.assertEqual('1', result['TEST2'])
        print result

        result = collapseAnnotation(self.ann_dict_1, 'max', fields='TEST2')
        self.assertEqual('2', result['TEST2'])
        print result

        result = collapseAnnotation(self.ann_dict_1, 'sum', fields='TEST2')
        self.assertEqual('3', result['TEST2'])
        print result
