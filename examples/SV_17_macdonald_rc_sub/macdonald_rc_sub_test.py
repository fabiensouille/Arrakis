#!/usr/bin/env python
import os
import unittest
from macdonald_vnv_1 import macdonald_1
from macdonald_vnv_2 import macdonald_2
from macdonald_vnv_3 import macdonald_3

class Checkmacdonald_rc_sub(unittest.TestCase):

    def test_macdonald_1(self):
        case = macdonald_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_macdonald_2(self):
        case = macdonald_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_macdonald_3(self):
        case = macdonald_3()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
