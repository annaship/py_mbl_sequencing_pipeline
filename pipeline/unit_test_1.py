import random
import unittest
#import sys
#from mock import Mock 
import db_upload as dbup

class MyConnectionTestCase(unittest.TestCase): 
    def test_execute_fetch_select(self): 
        mcn = dbup.MyConnection(host = "vampsdev", db = "test")
        table_name = "run_info_ill"
        id_name = table_name + "_id"
        sql = "select %s from %s where %s = 1" % (id_name, table_name, id_name)
        res = mcn.execute_fetch_select(sql)
        self.assertEqual(int(res[0][0]), 1)

if __name__ == '__main__':
    unittest.main()

