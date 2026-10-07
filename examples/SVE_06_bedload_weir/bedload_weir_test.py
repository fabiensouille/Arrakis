#!/usr/bin/env python
import os
import unittest
from bedload_weir_vnv_1 import bedload_weir_1

class Checkbedload_weir(unittest.TestCase):

    def test_bedload_weir_1(self):
        case = bedload_weir_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
