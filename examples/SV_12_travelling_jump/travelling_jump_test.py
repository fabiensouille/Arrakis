#!/usr/bin/env python
import os
import unittest
from travelling_jump_vnv_1 import travelling_jump_1
from travelling_jump_vnv_2 import travelling_jump_2

class Checktravelling_jump(unittest.TestCase):

    def test_travelling_jump_1(self):
        case = travelling_jump_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_travelling_jump_2(self):
        case = travelling_jump_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
