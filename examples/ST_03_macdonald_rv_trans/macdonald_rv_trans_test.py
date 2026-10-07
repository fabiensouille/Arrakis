#!/usr/bin/env python
import os
import unittest
from macdonald_vnv_1 import macdonald_1

class Checkmacdonald_rv_trans(unittest.TestCase):

    def test_macdonald_1(self):
        case = macdonald_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
