#!/usr/bin/env python
import os
import unittest
from vanrijn_vnv import vanrijn

class Checkbedload_dambreak(unittest.TestCase):

    def test_vanrijn_1(self):
        case = vanrijn()
        case.pre()
        case.run()
        case.check()
        case.post(show=False)

if __name__ == "__main__":
    unittest.main()
