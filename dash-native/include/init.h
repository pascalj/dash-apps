#pragma once

#include <libdash.h>
#include <array>

using md_type = double;

struct Float3D {
  md_type x;
  md_type y;
  md_type z;
};

struct Int3D {
  uint32_t x;
  uint32_t y;
  uint32_t z;
};

struct Input {
  Input(std::array<size_t, 3> size)
    : problem_size(size)
  {
  }

  std::array<size_t, 3> problem_size;
  float                 initial_temp = 1.44;
  float                 density      = 0.8442;
  float                 dt           = 0.005;
  float                 dtforce      = 0.0025;
  float                 force_cut    = 2.5;
  float                 cutneigh     = 2.8;
  float                 cutneighsq   = cutneigh * cutneigh;
  float                 neigh_every  = 20;
  float                 termo_every  = 100;
};

struct Config {
  using pattern_t = dash::BlockPattern<3>;

  // initialize config by generator
  Config() : input{{size, size, size}} {
    const float lattice = std::pow(4.0 / input.density, 1.0 / 3.0);
    for(uint8_t i = 0; i < 3; i++) {
      const auto prob_dim = input.problem_size[i];
      const auto lbox     = prob_dim * lattice;
      box[i]              = lbox;
      volume *= lbox;
      num_atoms *= prob_dim;
      num_bins[i] = 5.0 / 6.0 * prob_dim;
    }
    boxhi = box;
    // 4 atom types
    num_atoms *= 4;

  }

  uint32_t num_steps = 100;
  uint32_t size      = 32;
  uint32_t num_bins[3];

  // input variables
  Input input;

  std::array<md_type, 3> box, boxlo, boxhi;
  md_type     volume = 1.0;

  uint32_t num_atoms = 1;
  md_type  mass = 1.0;

  constexpr static md_type small  = 1.0e-6;
  constexpr static md_type factor = 0.999;
};

// RNG
const uint32_t IA   = 16807;
const uint32_t IM   = 2147483647;
const float    AM   = 1.0 / IM;
const uint32_t IQ   = 127773;
const uint32_t IR   = 2836;
const uint32_t MASK = 123459876;
