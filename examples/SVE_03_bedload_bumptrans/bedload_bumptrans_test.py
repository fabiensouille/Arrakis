#!/usr/bin/env python
import os
import unittest
from bedload_bumptrans_vnv_1 import bedload_bumptrans_1
from bedload_bumptrans_vnv_2 import bedload_bumptrans_2

class Checkbedload_bumptrans(unittest.TestCase):

    def test_bedload_bumptrans_1(self):
        case = bedload_bumptrans_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_bumptrans_2(self):
        case = bedload_bumptrans_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
