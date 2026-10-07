#!/usr/bin/env python
import os
import unittest
from tracer_multiple_vnv_1 import TracerMultiple_1

class CheckTracerMultiple(unittest.TestCase):

    def test_tracer_multiple_1(self):
        case = TracerMultiple_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)
        

if __name__ == "__main__":
    unittest.main()
