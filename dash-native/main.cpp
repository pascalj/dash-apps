#include <libdash.h>

#include <init.h>
#include <thermo.h>
#include <atom.h>
#include <force.h>

int main(int argc, char **argv) {
  dash::init(&argc, &argv);

  Config config;
  ForceLJ force(config.input.cutneigh);

  if (dash::myid() == 0) {
    std::cout << config;
  }

  auto termo = Thermo{config};
  Atoms atoms(config);

  if(dash::myid() == 0) {
    atoms.create_atoms();
  }

  // buildNeighbors();

  /* printSim(); */

  /* compute(true); */

  // timer start
  
  /* computeThermo(0, master); */


  // master timer start
 
  /* run(master); */

  // master timer end

  /* compute(true); */

  /* computeThermo(-1, master); */

  dash::finalize();
}
