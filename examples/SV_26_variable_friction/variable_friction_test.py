#!/usr/bin/env python
import os
import unittest
from variable_friction_vnv_1 import friction_1

class CheckVariableFriction(unittest.TestCase):

    def test_friction_1(self):
        case = friction_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
