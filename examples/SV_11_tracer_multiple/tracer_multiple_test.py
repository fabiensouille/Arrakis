#!/usr/bin/env python
import os
import unittest
from tracer_multiple_vnv_1 import TracerMultiple_1
from tracer_multiple_vnv_2 import TracerMultiple_2

class CheckTracerMultiple(unittest.TestCase):

    def test_tracer_multiple_1(self):
        case = TracerMultiple_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)
        
    def test_tracer_multiple_2(self):
        case = TracerMultiple_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
