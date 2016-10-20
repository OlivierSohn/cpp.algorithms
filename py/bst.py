
class BSTNode:
    """
    A node of a binary search tree

    """

    def __init__(self, key):
        self.left = None
        self.right = None
        self.key = key

    def traverse_range(self, traversal):
        left_traversed = False
        done = False

        if traversal.not_started():
            # look for the lower bound
            lower_bound = traversal.lower()
            left_traversed = True
            if lower_bound < self.key:
                if self.left is not None:
                    for n in self.left.traverse_range(traversal):
                        yield n
                if traversal.not_started():
                    traversal.start()
            elif lower_bound == self.key:
                traversal.start()
            else:
                if self.right is not None:
                    for n in self.right.traverse_range(traversal):
                        yield n
                done = True

        if (not done) and traversal.in_progress():
            if (self.left) and not left_traversed:
                for n in self.left.traverse_range(traversal):
                    yield n

            if traversal.in_progress():
                upper_bound = traversal.upper()
                if self.key > upper_bound:
                    traversal.stop()
                else:
                    yield self

                    if self.right is not None:
                        for n in self.right.traverse_range(traversal):
                            yield n

    def traverse(self):
        if self.left is not None:
            for n in self.left.traverse():
                yield n

        yield self

        if self.right is not None:
            for n in self.right.traverse():
                yield n

    def lower_bound(self, key):
        """
        returns the first element whose key is not considered to go before key

        """

        if key < self.key:
            if self.left:
                return self.left.lower_bound(key)
            return self
        if key > self.key:
            if self.right:
                return self.right.lower_bound(key)
            return
        return self

    def list_right(self, key_max):
        yield


    def insert(self, key):
        # we avoid recursion

        where = self
        while True:
            where = where._insert(key)
            if where is None:
                return

    def _insert(self, key):
        if key < self.key:
            if self.left is None:
                self.left = BSTNode(key)
            else:
                return self.left
        elif key == self.key:
            raise KeyError('duplicate keys are not supported in a BST')
        else:
            if self.right is None:
                self.right = BSTNode(key)
            else:
                return self.right

    def remove(self, key, parent):
        """
        Tries to removes a node by key.

        Returns True if the node was removed.
        Returns False if self should be removed by the caller

        """

        if key < self.key:
            self.left.remove(key, parent=self)
        elif key > self.key:
            self.right.remove(key, parent=self)
        else:
            self.remove_self(parent=parent)

    def remove_self(self, parent=None):
        """
        Removes self from the binary search tree structure

        :param parent: the parent of self or None if self is the root
        :return: the replacement for self
        """

        if (self.right is None) or (self.left is None):
            self_replacement = self.one_child()
            if parent is None:
                return self_replacement

            if parent.left == self:
                parent.left = self_replacement
            else:
                assert parent.right == self
                parent.right = self_replacement
        else:
            self.key = self.right.min_key()
            self.right.remove(self.key, parent=self)

        return self

    def one_child(self):
        return self.right if (self.left is None) else self.left

    def min_key(self):
        left_most = self
        while left_most.left is not None:
            left_most = left_most.left
        return left_most.key

    def count_nodes(self):
        n = 1
        if self.left:
            n += self.left.count_nodes()
        if self.right:
            n += self.right.count_nodes()
        return n

    def height(self):
        child_height = 0
        if self.left:
            child_height = self.left.count_nodes()
        if self.right:
            child_height = max(child_height, self.right.count_nodes())

        return 1 + child_height

_TRAVERSAL_NOT_STARTED = 0
_TRAVERSAL_IN_PROGRESS = 1
_TRAVERSAL_DONE = 2

class RangeTraversal(object):

    def __init__(self, min_key, max_key):
        self.min_key = min_key
        self.max_key = max_key
        self.state = _TRAVERSAL_NOT_STARTED

    def not_started(self):
        return self.state == _TRAVERSAL_NOT_STARTED

    def in_progress(self):
        return self.state == _TRAVERSAL_IN_PROGRESS

    def stop(self):
        assert self.state != _TRAVERSAL_DONE
        self.state = _TRAVERSAL_DONE

    def start(self):
        assert self.state == _TRAVERSAL_NOT_STARTED
        self.state = _TRAVERSAL_IN_PROGRESS

    def lower(self):
        return self.min_key

    def upper(self):
        return self.max_key


class BST(object):
    """ A binary search tree

    """

    def __init__(self):
        self.root = None

    def append(self, key):
        """ Appends a key

        """

        if self.root is None:
            self.root = BSTNode(key)
        else:
            self.root.insert(key)

    def remove(self, key):
        """ Removes a key

        """

        assert self.root is not None, 'cannot remove from an empty tree'

        if self.root.key == key:
            self.root = self.root.remove_self()
        else:
            self.root.remove(key, parent=None)

    def traverse_range(self, key_min, key_max):
        """
        Generator for a range of elements

        """

        assert key_min <= key_max, 'order of keys is not consistent'

        if self.root is not None:
            t = RangeTraversal(key_min, key_max)
            for n in self.root.traverse_range(t):
                yield n.key

    def fill_rate(self):
        return float(len(self)) / float(self.theoretical_length_given_height())

    def theoretical_length_given_height(self):
        height = self.height()
        if height == 0:
            return 0
        return 1 + pow(2, height - 1)

    def height(self):
        if self.root is None:
            return 0
        return self.root.height()

    def __len__(self):
        if self.root is None:
            return 0
        return self.root.count_nodes()
