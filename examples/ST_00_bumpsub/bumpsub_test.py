#!/usr/bin/env python
import os
import unittest
from bumpsub_vnv_1 import bumpsub_1

class CheckBumpsub(unittest.TestCase):

    def test_bumpsub_1(self):
        case = bumpsub_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
