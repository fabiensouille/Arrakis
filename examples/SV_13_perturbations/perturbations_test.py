#!/usr/bin/env python
import os
import unittest
from perturbations_vnv_1 import perturbations_1
from perturbations_vnv_2 import perturbations_2

class Checkperturbations(unittest.TestCase):

    def test_perturbations_1(self):
        case = perturbations_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_perturbations_2(self):
        case = perturbations_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
