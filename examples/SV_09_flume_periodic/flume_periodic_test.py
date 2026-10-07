#!/usr/bin/env python
import os
import unittest
from flume_vnv_1 import flume_1
from flume_vnv_2 import flume_2

class Checkflume(unittest.TestCase):

    def test_flume_1(self):
        case = flume_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

    def test_flume_2(self):
        case = flume_2()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)


if __name__ == "__main__":
    unittest.main()
