#!/usr/bin/env python
import os
import unittest
from solitary_wave_vnv_1 import solitary_wave_1
from solitary_wave_vnv_2 import solitary_wave_2

class Checksolitary_wave(unittest.TestCase):

    def test_solitary_wave_1(self):
        case = solitary_wave_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_solitary_wave_2(self):
        case = solitary_wave_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
