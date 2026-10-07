#!/usr/bin/env python
import os
import unittest
from friction_vnv_1 import friction_1
from friction_vnv_2 import friction_2
from friction_vnv_3 import friction_3
from friction_vnv_4 import friction_4
from friction_vnv_5 import friction_5

class Checkfriction(unittest.TestCase):

    def test_friction_1(self):
        case = friction_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_friction_2(self):
        case = friction_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)
        
    def test_friction_3(self):
        case = friction_3()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_friction_4(self):
        case = friction_4()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_friction_5(self):
        case = friction_5()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
