#!/usr/bin/env python
import os
import unittest
from bedload_soni_vnv_1 import bedload_soni_1

class Checkbedload_soni(unittest.TestCase):

    def test_bedload_soni_1(self):
        case = bedload_soni_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
