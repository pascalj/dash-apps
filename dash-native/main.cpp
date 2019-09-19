#include <libdash.h>

#include <init.h>

int main(int argc, char **argv) {
  dash::init(&argc, &argv);

  Config config;

  // initThermo();

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
