#!/usr/bin/env python
import os
import unittest
from deposition_vnv import deposition

class Checkbedload_dambreak(unittest.TestCase):

    def test_deposition_1(self):
        case = deposition()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
