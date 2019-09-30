#pragma once

#include "init.h"
#include "force.h"

struct Thermo {
  Thermo() = delete;

  Thermo(Config& config)
    : config(config)
    , mvv2e(1.0)
    , dof_boltz(config.num_atoms)
    , temp_scale(mvv2e / dof_boltz)
    , press_scale(1.0 / 3 / config.volume)
    , eng_scale(0.5)
  {
    config.input.temp_scale = temp_scale;
  }

  Config &config;
  float mvv2e = 1.0;
  float dof_boltz;
  float temp_scale;
  float press_scale;
  float eng_scale;

  float temperature() const {
    // TODO: implement
    return 0;
  }

  void compute(Force &force, int step) {
    if (step > 0 && step % config.input.thermo_every != 0) {
      return;
    }

    auto temp  = temperature();
    auto eng   = force.eng_vdwl * eng_scale / config.num_atoms;
    auto press = (temp * dof_boltz + force.virial) * press_scale;

    
  }

};
