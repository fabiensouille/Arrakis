#!/usr/bin/env python
import os
import unittest
from bedload_antidune_vnv_1 import bedload_antidune_1
from bedload_antidune_vnv_2 import bedload_antidune_2

class Checkbedload_antidune(unittest.TestCase):

    def test_bedload_antidune_1(self):
        case = bedload_antidune_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_antidune_2(self):
        case = bedload_antidune_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
