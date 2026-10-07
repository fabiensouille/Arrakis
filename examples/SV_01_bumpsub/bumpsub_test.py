#!/usr/bin/env python
import os
import unittest
from bumpsub_vnv_1 import bumpsub_1
from bumpsub_vnv_2 import bumpsub_2

class CheckBumpsub(unittest.TestCase):

    def test_bumpsub_1(self):
        case = bumpsub_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bumpsub_2(self):
        case = bumpsub_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
