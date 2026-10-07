#!/usr/bin/env python
import os
import unittest
from run_up_vnv_1 import run_up_1
from run_up_vnv_2 import run_up_2

class Checkrun_up(unittest.TestCase):

    def test_run_up_1(self):
        case = run_up_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_run_up_2(self):
        case = run_up_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
