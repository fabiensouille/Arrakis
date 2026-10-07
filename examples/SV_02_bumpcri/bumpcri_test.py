#!/usr/bin/env python
import os
import unittest
from bumpcri_vnv_1 import bumpcri_1
from bumpcri_vnv_2 import bumpcri_2

class CheckBumpcri(unittest.TestCase):

    def test_bumpcri_1(self):
        case = bumpcri_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bumpcri_2(self):
        case = bumpcri_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
