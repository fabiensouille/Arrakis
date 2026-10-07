#!/usr/bin/env python
import os
import unittest
from tracer_canal_vnv import TracerCanal

class CheckTracerCanal(unittest.TestCase):

    def test_tracer_canal_1(self):
        case = TracerCanal()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
