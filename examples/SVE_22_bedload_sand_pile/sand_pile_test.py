#!/usr/bin/env python
import os
import unittest
from sand_pile_vnv_1 import sand_pile_1

class Checksandpile(unittest.TestCase):

    def test_sand_pile_1(self):
        case = sand_pile_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
