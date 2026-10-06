#!/usr/bin/env python
import os
import unittest
from bedload_dambreak_vnv_1 import bedload_dambreak_1
from bedload_dambreak_vnv_2 import bedload_dambreak_2

class Checkbedload_dambreak(unittest.TestCase):

    def test_bedload_dambreak_1(self):
        case = bedload_dambreak_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_dambreak_2(self):
        case = bedload_dambreak_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
