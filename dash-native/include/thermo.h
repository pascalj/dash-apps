#pragma once

#include "init.h"
#include "force.h"

struct Stats {
  Stats(md_type temp, md_type eng, md_type press)
    : temp(temp)
    , eng(eng)
    , press(press)
  {
  }

  md_type temp;
  md_type eng;
  md_type press;
};

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
  md_type mvv2e = 1.0;
  md_type dof_boltz;
  md_type temp_scale;
  md_type press_scale;
  md_type eng_scale;

  double temperature() const {
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

    stats.emplace_back(temp, eng, press);
  }

  std::vector<Stats> stats;
};
