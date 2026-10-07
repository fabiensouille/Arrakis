#!/usr/bin/env python
import os
import unittest
from bedload_stream_vnv_1 import bedload_stream_1
from bedload_stream_vnv_2 import bedload_stream_2

class Checkbedload_stream(unittest.TestCase):

    def test_bedload_stream_1(self):
        case = bedload_stream_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_bedload_stream_2(self):
        case = bedload_stream_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
