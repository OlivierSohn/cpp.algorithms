
class MinHeap(object):
    """ List-based min heap

    """

    def __init__(self, elements=None):
        if elements is None:
            elements = []
        self.array = elements

        if len(self) >= 2:
            self._sort()

    def append(self, value):
        """ Appends a value and maintain the min heap property.

            The value is inserted as a leaf, and "bubbled up" by successive swaps
            with its immediate parent.

        """
        self.array.append(value)

        my_index = len(self.array)
        while True:
            if my_index == 1:
                # 1 is the index of the root : we can't bubble up any more
                return
            parent_index = my_index // 2
            parent_value = self.array[parent_index - 1]
            if value >= parent_value:
                return
            self._swap(parent_index-1, my_index-1)
            my_index = parent_index

    def empty(self):
        return self.size() == 0

    def size(self):
        return len(self)

    def __len__(self):
        return len(self.array)

    def min(self):
        return self.array[0]

    def pop(self):
        # the min element is at position 0

        # swap first and last elements
        self._swap(0,
                   -1)

        # pop the last element (which is now the min element)
        elem = self.array.pop()

        # sort the first element
        self._min_heapify(0)

        return elem

    def _swap(self, index1, index2):
        self.array[index1], self.array[index2] = self.array[index2], self.array[index1]

    def _min_heapify(self, index_root):
        index_root += 1

        left_index = 2 * index_root
        size = self.size()
        if left_index > size:
            # no left child
            return

        left_value = self.array[ left_index - 1 ]
        root_value = self.array[ index_root - 1 ]
        right_index = left_index + 1
        if right_index > size:
            # no right child
            if left_value >= root_value:
                return

            self._swap(left_index - 1,
                       index_root - 1)
            self._min_heapify(left_index - 1)
            return
        # right child
        right_index = left_index + 1
        right_value = self.array[ right_index - 1 ]

        if left_value < right_value:
            if left_value >= root_value:
                return

            self._swap(left_index - 1,
                       index_root - 1)
            self._min_heapify(left_index - 1)
            return

        if right_value >= root_value:
            return
        self._swap(right_index - 1,
                   index_root - 1)

        self._min_heapify(right_index - 1)


    def _sort(self):
        index = (len(self) - 1) // 2
        while index != -1:
            self._min_heapify(index)
            index -= 1
