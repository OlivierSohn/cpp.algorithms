import unittest

from heap import MinHeap

class TestHeap(unittest.TestCase):

    def test_empty(self):
        mh = MinHeap()
        self.assertTrue(mh.empty())

        mh.append(1)
        self.assertFalse(mh.empty())

        mh.pop()
        self.assertTrue(mh.empty())

    def test_len(self):
        mh = MinHeap()
        self.assertTrue(len(mh) == 0)

        mh.append(1)
        self.assertTrue(len(mh) == 1)

        mh.append(1)
        self.assertTrue(len(mh) == 2)

        mh.pop()
        self.assertTrue(len(mh) == 1)

        mh.pop()
        self.assertTrue(len(mh) == 0)

    def test_append(self):
        mh = MinHeap()
        self.assertEqual(0, mh.size())

        mh.append(1)
        self.assertEqual(1, mh.size())

        mh.append(1)
        self.assertEqual(2, mh.size())

        mh.append(2)
        self.assertEqual(3, mh.size())

    def test_min(self):
        mh = MinHeap()
        with self.assertRaises(IndexError):
            mh.min()

        mh.append(1)
        self.assertEqual(1, mh.min())
        self.assertEqual(1, mh.min())

        mh.append(1)
        self.assertEqual(1, mh.min())

        mh.pop()
        self.assertEqual(1, mh.min())

        mh.pop()
        with self.assertRaises(IndexError):
            mh.min()

        mh.append(3)
        self.assertEqual(3, mh.min())

        mh.append(1)
        self.assertEqual(1, mh.min())

        mh.pop()
        self.assertEqual(3, mh.min())

    def test_pop(self):
        mh = MinHeap()
        with self.assertRaises(IndexError):
            mh.pop()

        mh.append(1)
        self.assertEqual(1, mh.pop())
        with self.assertRaises(IndexError):
            mh.pop()

        mh.append(9)
        mh.append(6)
        mh.append(5)
        mh.append(3)
        self.assertEqual(3, mh.pop())
        self.assertEqual(5, mh.pop())
        self.assertEqual(6, mh.pop())
        self.assertEqual(9, mh.pop())
        with self.assertRaises(IndexError):
            mh.pop()

if __name__ == '__main__':
    unittest.main()
