#!/usr/bin/env python
import os
import unittest
from dambreak_constriction_vnv_1 import constriction_1

class CheckDambreak(unittest.TestCase):

    def test_constriction_1(self):
        case = constriction_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
