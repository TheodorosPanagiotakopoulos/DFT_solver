#pragma once
#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>

namespace dft { namespace units {

struct au {
  static constexpr double hartree_to_ev = 27.211386245988;
  static constexpr double bohr_to_ang = 0.529177210903;
};

inline double ang_to_bohr(double ang) { return ang / au::bohr_to_ang; }
inline double bohr_to_ang(double bohr){ return bohr * au::bohr_to_ang; }

}}
