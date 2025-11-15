#include <iostream>
#include "dft/grid.hpp"
#include "dft/potential.hpp"
#include "dft/scf.hpp"
#include "dft/io.hpp"

using namespace dft;

int main(){
  Grid3D g({20,20,20},{10.0,10.0,10.0});
  pot::HarmonicOscillator vext(0.5);
  SCFOptions opt;
  opt.nelectrons=2;
  opt.max_iters=25;

  auto res = run_scf(g,vext,opt);
  std::cout<<"\nConverged in "<<res.iters<<" iterations\n";
  std::cout<<"Total Energy = "<<res.Etot<<"\n";
  io::write_scalar("density.txt",res.density,g.n[0],g.n[1],g.n[2]);
}
