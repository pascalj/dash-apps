#pragma once

#include <libdash.h>
#include "init.h"

using Velocity = Float3D;
using Force    = Float3D;

struct Atom {
  Velocity v;
  Force    f;
};

struct Atoms {
  Atoms() = delete;

  Atoms(const Config& cfg)
  {
    pos.allocate(cfg.pos_spec());
    space.allocate(cfg.space_spec());
    bins.allocate(cfg.bin_spec());
  }

  dash::NArray<Float3D, 4> pos;
  dash::NArray<Float3D, 4> space;
  // actual storage of atoms in bins
  dash::NArray<Atom, 4>    bins;
};
