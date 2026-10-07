#!/usr/bin/env python
import os
import unittest
from dambreak_triangular_obstacle_vnv_1 import dambreak_triangular_obstacle_1
from dambreak_triangular_obstacle_vnv_2 import dambreak_triangular_obstacle_2

class Checkdambreak_triangular_obstacle(unittest.TestCase):

    def test_dambreak_triangular_obstacle_1(self):
        case = dambreak_triangular_obstacle_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_dambreak_triangular_obstacle_2(self):
        case = dambreak_triangular_obstacle_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
