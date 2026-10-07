#!/usr/bin/env python
import os
import unittest
from slope_effect_vnv import slope_effect

class Checkbedload_dambreak(unittest.TestCase):

    def test_slope_effect_1(self):
        case = slope_effect()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
