#pragma once

#include <init.h>

struct Neighbor {
  int x;
  int y;
  int z;

  int index;
};

struct Neighbors {
  Neighbors() = delete;

  Neighbors(const Config& cfg)
  {
  }

  dash::NArray<Neighbor, 2> neighs;
  uint32_t                  ncount = 0;
};
