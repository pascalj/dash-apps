#pragma once

#include "init.h"

struct Thermo {
  Thermo() = delete;

  Thermo(Config& config)
    : mvv2e(1.0)
    , dof_boltz(config.num_atoms)
    , temp_scale(mvv2e / dof_boltz)
    , press_scale(1.0 / 3 / config.volume)
    , eng_scale(0.5)
  {
  }

  float mvv2e = 1.0;
  float dof_boltz;
  float temp_scale;
  float press_scale;
  float eng_scale;
};
