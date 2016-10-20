import unittest
from random import shuffle, seed

from bst import BST


class TestBst(unittest.TestCase):

    def test_one_element(self):
        b = BST()
        b.append(1)
        with self.assertRaises(KeyError):
            b.append(1)

        b.remove(1)
        with self.assertRaises(AssertionError):
            b.remove(1)


    def test_multiple_elements(self):
        b = BST()

        for i in xrange(0,10):
            b.append(i)

        self.assertIsNotNone(b.root)

        for i in xrange(0,10):
            b.remove(i)

        self.assertIsNone(b.root)


    def test_range_simple(self):
        b = BST()

        l = list(b.traverse_range(0, 0))
        self.assertEqual(0, len(l))

        b.append(0)

        l = list(b.traverse_range(0, 0))
        self.assertEqual(1, len(l))

        b.append(1)
        b.append(-1)

        l = list(b.traverse_range(0, 0))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(-0.1, 0.1))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(-0.1, 0))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(0, 0.1))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(-1, -1))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(-1.1, -0.9))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(-1.1, -1))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(-1, -0.9))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(1, 1))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(0.9, 1.1))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(1, 1.1))
        self.assertEqual(1, len(l))
        l = list(b.traverse_range(0.9, 1))
        self.assertEqual(1, len(l))

        l = list(b.traverse_range(0, 1))
        self.assertEqual(2, len(l))
        l = list(b.traverse_range(-0.1, 1))
        self.assertEqual(2, len(l))
        l = list(b.traverse_range(0, 1.1))
        self.assertEqual(2, len(l))
        l = list(b.traverse_range(-0.1, 1.1))
        self.assertEqual(2, len(l))

        l = list(b.traverse_range(-1, 1))
        self.assertEqual(3, len(l))
        l = list(b.traverse_range(-1.1, 1))
        self.assertEqual(3, len(l))
        l = list(b.traverse_range(-1.1, 1.1))
        self.assertEqual(3, len(l))
        l = list(b.traverse_range(-1, 1.1))
        self.assertEqual(3, len(l))

        l = list(b.traverse_range(-10.6, 1.2))
        self.assertEqual(3, len(l))

    def test_range_complex(self):

        seed(0)
        for _ in range(10):
            keys = [i for i in range(1000)]
            shuffle(keys)

            self.with_keys(keys)

    def with_keys(self, keys):
        length_keys = len(keys)

        b = BST()
        for key in keys:
            b.append(key)

        result_list = list(b.traverse_range(0, 0))
        self.verify_length(result_list, length_keys, 0, 0)

        for min_ in [-10, 0, 1, 2, 10, 100, 500, 998, 999, 1000, 1001, 1010]:
            for max_ in [min_, min_+1, min_+2, min_+3, min_+100, min_+1000]:

                result_list = list(b.traverse_range(min_, max_))
                self.verify_length(result_list, length_keys, min_, max_)
                prev = None
                for i in result_list:
                    if prev is None:
                        prev = i
                        continue
                    self.assertEqual(prev, i-1)
                    prev = i

    def verify_length(self, result_list, length_keys, min_, max_):
        if min_ < 0:
            min_ = 0
        elif min_ >= length_keys:
            min_ = length_keys

        if max_ >= length_keys:
            max_ = length_keys - 1

        expected_length = max_ - min_ + 1
        if expected_length < 0:
            expected_length = 0

        self.assertEqual(expected_length, len(result_list))


if __name__ == '__main__':
    unittest.main()

