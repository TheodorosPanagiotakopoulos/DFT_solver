#pragma once
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>

namespace dft { namespace io {

inline void write_scalar(const std::string& file,
                         const std::vector<double>& data,
                         size_t nx,size_t ny,size_t nz)
{
  std::ofstream ofs(file);
  ofs << std::setprecision(10);
  for(size_t i=0;i<data.size();++i) ofs << data[i] << "\n";
}

}
}
