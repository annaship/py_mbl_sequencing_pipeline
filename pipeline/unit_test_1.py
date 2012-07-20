import random
import unittest
#import sys
#from mock import Mock 
import db_upload as dbup

class TestSequenceFunctions(unittest.TestCase):

    def setUp(self):
        self.seq = list(range(10))

    def test_shuffle(self):
        # make sure the shuffled sequence does not lose any elements
        random.shuffle(self.seq)
        self.seq.sort()
        self.assertEqual(self.seq, list(range(10)))

        # should raise an exception for an immutable sequence
        self.assertRaises(TypeError, random.shuffle, (1,2,3))

    def test_choice(self):
        element = random.choice(self.seq)
        self.assertTrue(element in self.seq)

    def test_sample(self):
        with self.assertRaises(ValueError):
            random.sample(self.seq, 20)
        for element in random.sample(self.seq, 5):
            self.assertTrue(element in self.seq)
            
class MyConnectionTestCase(unittest.TestCase): 
    def test_execute_fetch_select(self): 
        #set up the mock objects 
#        mockCursor = Mock()
#        mockDB = Mock( { "cursor" : mockCursor } ) 
        #call the function to be tested
        # persistData(testData, mockDB) 
        #test the correct calls were made on the database 
        # print dir(mockDB)
        # ['__class__', '__cmp__', '__contains__', '__delattr__', '__delitem__', '__doc__', '__eq__', '__format__', '__ge__', 
        # '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__iter__', '__le__', '__len__', '__lt__', '__ne__',
        # '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setitem__', '__sizeof__', '__str__', '__subclasshook__', 
        # 'assert_any_call', 'assert_called_once_with', 'assert_called_with', 'assert_has_calls', 'attach_mock', 'call_args', 'call_args_list', 
        # 'call_count', 'called', 'clear', 'configure_mock', 'copy', 'fromkeys', 'get', 'has_key', 'items', 'iteritems', 'iterkeys', 'itervalues', 
        # 'keys', 'method_calls', 'mock_add_spec', 'mock_calls', 'pop', 'popitem', 'reset_mock', 'return_value', 'setdefault', 'side_effect', 'update',
        # 'values', 'viewitems', 'viewkeys', 'viewvalues']
        mcn = dbup.MyConnection()
        sql = "select run_info_ill_id from run_info_ill where run_info_ill_id = 1"
        res = mcn.execute_fetch_select(sql)
        self.assertEqual(int(res[0][0]), 1)

        # objects mockDB.mockCheckCall(0, 'cursor')
        # mockCursor.mockCheckCall(0, 
        #                          'execute', 
        #                          '...some SQL...',
        #                           someData) 
        # mockCursor.mockCheckCall(1, 
        #                          'execute', 
        #                          '...more SQL...', 
        #                           moreData) 
        # mockDB.mockCheckCall(1, 'commit')
        # mockDB.mockCheckCall(2, 'close')


if __name__ == '__main__':
    unittest.main()

