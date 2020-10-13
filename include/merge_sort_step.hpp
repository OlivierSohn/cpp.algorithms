namespace imajuscule {

// Performs one "merge" step of the merge step sort.
//
// Preconditions:
// - ranges [it1, end1) and [it2, end2) are sorted according to 'f'
template<typename It, typename F, typename Container>
void merge_sort_step(It it1,
                     It const &end1,
                     It it2,
                     It const &end2,
                     F && f,
                     Container & res) {
  res.clear();
  res.reserve(std::distance(it1, end1) + std::distance(it2, end2));
  
  if (it1 == end1) {
    for (; it2 != end2; ++it2) {
      res.push_back(*it2);
    }
    return;
  }
  if (it2 == end2) {
    for (; it1 != end1; ++it1) {
      res.push_back(*it1);
    }
    return;
  }

  while (true) {
    if (f(*it1) < f(*it2)) {
      res.push_back(*it1);
      ++it1;
      if (it1 == end1) {
        for (; it2 != end2; ++it2) {
          res.push_back(*it2);
        }
        return;
      }
    } else {
      res.push_back(*it2);
      ++it2;
      if (it2 == end2) {
        for (; it1 != end1; ++it1) {
          res.push_back(*it1);
        }
        return;
      }
    }
  }
}

}
