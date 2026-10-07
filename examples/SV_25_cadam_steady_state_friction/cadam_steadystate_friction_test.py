#!/usr/bin/env python
import os
import unittest
from cadam_steadystate_vnv_0 import steadystate_0
from cadam_steadystate_vnv_1 import steadystate_1

class CheckSteadyState(unittest.TestCase):

    def test_steadystate_0(self):
        case = steadystate_0()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_steadystate_1(self):
        case = steadystate_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
