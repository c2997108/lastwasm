// Author: Martin C. Frith 2019
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef MCF_CONTIGUOUS_QUEUE
#define MCF_CONTIGUOUS_QUEUE

#include <vector>
#include <stddef.h>  // size_t

namespace mcf {

template <typename T> class ContiguousQueue {
public:
  void clear() {
    v.clear();
    start = 0;
  }

  size_t size() const { return start + v.size(); }

  void push(const T &item, size_t newStart) {
    size_t oldSize = newStart - start;
    if (oldSize > v.size() / 2) {
      v.erase(v.begin(), v.begin() + oldSize);
      start = newStart;
    }
    v.push_back(item);
  }

  const T &fromEnd(int n) const { return v.end()[-n]; }

private:
  std::vector<T> v;
  size_t start;
};

}

#endif
