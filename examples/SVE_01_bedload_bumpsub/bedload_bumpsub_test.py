#!/usr/bin/env python
import os
import unittest
from bedload_bumpsub_vnv_1 import bedload_bumpsub_1
from bedload_bumpsub_vnv_2 import bedload_bumpsub_2

class Checkbedload_bumpsub(unittest.TestCase):

    def test_bedload_bumpsub_1(self):
        case = bedload_bumpsub_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_bumpsub_2(self):
        case = bedload_bumpsub_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
