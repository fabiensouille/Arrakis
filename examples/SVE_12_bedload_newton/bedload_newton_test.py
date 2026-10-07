#!/usr/bin/env python
import os
import unittest
from bedload_newton_vnv_1 import bedload_newton_1

class Checkbedload_newton(unittest.TestCase):

    def test_bedload_newton_1(self):
        case = bedload_newton_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
