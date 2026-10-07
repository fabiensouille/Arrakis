#!/usr/bin/env python
import os
import unittest
from bedload_bumpcri_vnv_1 import bedload_bumpcri_1
from bedload_bumpcri_vnv_2 import bedload_bumpcri_2

class Checkbedload_bumpcri(unittest.TestCase):

    def test_bedload_bumpcri_1(self):
        case = bedload_bumpcri_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_bumpcri_2(self):
        case = bedload_bumpcri_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
