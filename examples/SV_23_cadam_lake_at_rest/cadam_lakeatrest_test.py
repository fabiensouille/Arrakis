#!/usr/bin/env python
import os
import unittest
from cadam_lakeatrest_vnv_0 import lakeatrest_0
from cadam_lakeatrest_vnv_1 import lakeatrest_1
from cadam_lakeatrest_vnv_2 import lakeatrest_2

class Checklakeatrest(unittest.TestCase):

    def test_lakeatrest_0(self):
        case = lakeatrest_0()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_lakeatrest_1(self):
        case = lakeatrest_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_lakeatrest_2(self):
        case = lakeatrest_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
