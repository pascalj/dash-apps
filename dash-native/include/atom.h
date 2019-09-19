#pragma once

#include "init.h"

using Velocity = Float3D;
using Force    = Float3D;

struct Atoms {
  Velocity v;
  Force f;
};
