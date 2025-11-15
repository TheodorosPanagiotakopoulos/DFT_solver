#pragma once
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace dft { namespace par {

template<typename F>
void for_index(size_t N, F f){
  tbb::parallel_for(tbb::blocked_range<size_t>(0,N),
    [&](const tbb::blocked_range<size_t>& r){
      for(size_t i=r.begin(); i<r.end(); ++i) f(i);
    });
}

}}
