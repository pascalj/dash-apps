#include <libdash.h>

#include <init.h>
#include <thermo.h>
#include <atom.h>
#include <neighbors.h>
#include <force.h>
#include <timer.h>
#include <mpi.h>

#include <algorithm>

int main(int argc, char **argv) {
  dash::init(&argc, &argv);

  Config  config;
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

  force.compute(atoms, neighbors, true);

  const size_t layer = 1;
      for (size_t i = 0; i < 6; i++) {
        for (size_t j = 0; j < 6; j++) {
          /* std::cout << std::setw(2) << atoms.per_bin[i][j][layer].get() << " "; */
          Float3D pos{0.0, 0.0, 0.0};
          for(size_t k = 0; k < atoms.per_bin[i][j][layer]; k++) {
            pos[0] += atoms.atoms[i][j][layer][k].get().pos[0];
            pos[1] += atoms.atoms[i][j][layer][k].get().pos[1];
            pos[2] += atoms.atoms[i][j][layer][k].get().pos[2];
            /* std::cout << atoms.atoms[i][j][0][k].get() << std::endl; */
          }
          auto per_bin = atoms.per_bin[i][j][layer].get();
          if(per_bin > 0) {
            printf("(%05.2f, %05.2f, %05.2f) ", pos[0] / per_bin, pos[1] / per_bin, pos[2] / per_bin);
          } else {
            printf("(...................) ");
          }
        }

        std::cout << "  |  ";
        for (size_t j = 0; j < config.num_bins[1]; j++) {
          std::cout << std::setw(2) << neighbors.neighs[i * config.num_bins[0] + j].size() << " ";
        }
        std::cout << std::endl;
      }
  // timer start
  thermo.compute(force, 0);

  // master timer start

  // main loop
  double start = MPI_Wtime();
  for(uint32_t step = 0; step < config.num_steps; step++) {

    /* t(start, "Iteration start"); */
    atoms.initial_integrate();
    /* t(start, "initial integrate done"); */

    if (step % config.input.neigh_every == 0) {
      neighbors.rebuild(atoms);
    } else {
      neighbors.update_positions(atoms);
    }

    /* t(start, "neighbors done"); */
    force.compute(atoms, neighbors, step % config.input.thermo_every == 0);

    /* t(start, "force compute done"); */
    atoms.final_integrate();

    /* t(start, "final integrate done"); */
    if(config.input.thermo_every > 0) {
      thermo.compute(force, step);
      t(start, "thermo compute done");
    }

    std::cout << "t " << atoms.temperature() << std::endl;
  }
  t(start, "end");

  dash::finalize();
}

