#!/usr/bin/env python
import os
import unittest
from bumptrans_vnv_1 import bumptrans_1
from bumptrans_vnv_2 import bumptrans_2

class Checkbumptrans(unittest.TestCase):

    def test_bumptrans_1(self):
        case = bumptrans_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bumptrans_2(self):
        case = bumptrans_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
