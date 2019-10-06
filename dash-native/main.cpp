#include <libdash.h>

#include <init.h>
#include <thermo.h>
#include <atom.h>
#include <neighbors.h>
#include <force.h>
#include <timer.h>
#include <mpi.h>

int main(int argc, char **argv) {
  dash::init(&argc, &argv);

  Config config;
  ForceLJ force(config.input.cutneigh);

  if (dash::myid() == 0) {
    std::cout << config;
  }

  auto thermo = Thermo{config};
  Atoms atoms(config);
  Neighbors neighbors(config);

  if(dash::myid() == 0) {
    atoms.create_atoms();
  }
  dash::barrier();
  atoms.create_velocity();

  neighbors.rebuild(atoms);

  /* printSim(); */

  force.compute(atoms, neighbors, true);

  // timer start
  
  thermo.compute(force, 0);

  /* int i = 0; */
  /* for(auto atom : atoms.atoms) { */
  /*   std::cout << "[" << i++ << "] " << (atom).get() <<std::endl; */
  /* } */
  /* return 1; */


  // master timer start

  // main loop
  double start = MPI_Wtime();
  for(uint32_t step = 0; step < config.num_steps; step++) {
    
    t(start, "Iteration start");
    atoms.initial_integrate();
    t(start, "initial integrate done");

    if (step % config.input.neigh_every == 0) {
      neighbors.rebuild(atoms);
    } else {
      neighbors.update_positions(atoms);
    }

    t(start, "neighbors done");
    force.compute(atoms, neighbors, step % config.input.thermo_every == 0);

    t(start, "force compute done");
    atoms.final_integrate();

    t(start, "final integrate done");
    if(config.input.thermo_every > 0) {
      thermo.compute(force, step);
      t(start, "thermo compute done");

    }

  }

  dash::finalize();
}
