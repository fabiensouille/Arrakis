#!/usr/bin/env python
import os
import unittest
from saint_matthieu_std import stmatthieu_1

class CheckStMatthieu(unittest.TestCase):

    def test_stmatthieu_1(self):
        case = stmatthieu_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
