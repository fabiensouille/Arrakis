#!/usr/bin/env python
import os
import unittest
from flume_vnv import flume

class Checkflume(unittest.TestCase):

    def test_dambreak(self):
        case = flume()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
