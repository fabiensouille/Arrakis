#!/usr/bin/env python
import os
import unittest
from channel_constriction_vnv_1 import channel_constriction_1

class Check_channel_constriction(unittest.TestCase):

    def test_channel_constriction_1(self):
        case = channel_constriction_1()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
