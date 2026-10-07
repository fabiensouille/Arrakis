#!/usr/bin/env python
import os
import unittest
from waves_over_a_bar_vnv_1 import waves_over_a_bar_1
from waves_over_a_bar_vnv_2 import waves_over_a_bar_2

class Checkwaves_over_a_bar(unittest.TestCase):

    def test_waves_over_a_bar_1(self):
        case = waves_over_a_bar_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_waves_over_a_bar_2(self):
        case = waves_over_a_bar_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
