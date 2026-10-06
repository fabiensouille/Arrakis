#!/usr/bin/env python
import os
import unittest
from dambreak_vnv_1 import Dambreak_1
from dambreak_vnv_2 import Dambreak_2
from dambreak_vnv_3 import Dambreak_3
from dambreak_vnv_4 import Dambreak_4

class CheckDambreak(unittest.TestCase):

    def test_dambreak_1(self):
        case = Dambreak_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_dambreak_2(self):
        case = Dambreak_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_dambreak_3(self):
        case = Dambreak_3()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
