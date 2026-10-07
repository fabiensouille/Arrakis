#!/usr/bin/env python
import os
import unittest
from thacker_vnv_1 import thacker_1
from thacker_vnv_2 import thacker_2

class Checkthacker(unittest.TestCase):

    def test_thacker_1(self):
        case = thacker_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_thacker_2(self):
        case = thacker_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
