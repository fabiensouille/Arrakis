#!/usr/bin/env python
import os
import unittest
from bedload_lakeatrest_vnv_1 import bedload_lakeatrest_1
from bedload_lakeatrest_vnv_2 import bedload_lakeatrest_2

class CheckbedloadLakeAtRest(unittest.TestCase):

    def test_bedload_lakeatrest_1(self):
        case = bedload_lakeatrest_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_lakeatrest_2(self):
        case = bedload_lakeatrest_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
