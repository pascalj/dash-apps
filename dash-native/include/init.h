#pragma once

#include <libdash.h>
#include <valarray>

using md_type = double;

struct Float3D {
  explicit Float3D(std::valarray<md_type> i)
    : x(i[0])
    , y(i[1])
    , z(i[2])
  {
  }

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
    , lattice(std::pow(4.0 / density, 1.0 / 3.0))
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
  int                   neigh_every  = 20;
  int                   thermo_every = 100;
  float                 lattice;
};

struct Config {
  using pattern_t = dash::BlockPattern<3>;

  // initialize config by generator
  Config()
    : input{{size, size, size}}
    , box(3)
    , boxlo(3)
    , boxhi(3)
  {
    const float lattice = input.lattice;
    for (uint8_t i = 0; i < 3; i++) {
      const auto prob_dim = input.problem_size[i];
      const auto lbox     = prob_dim * lattice;
      box[i]              = lbox;
      volume *= lbox;
      num_atoms *= prob_dim;
      num_bins[i]   = 5.0 / 6.0 * prob_dim;
      bin_size[i]   = box[i] / num_bins[i];
      bin_inv[i]    = 1.0 / bin_size[i];
      bin_needed[i] = input.cutneigh * bin_inv[i];
      if (bin_needed[i] * bin_size[i] < factor * input.cutneigh) {
        bin_needed[i] += 1;
      }
    }
    boxhi = box;
  }

  auto pos_spec() const
  {
    return dash::SizeSpec<4>(
        num_bins[0] + bin_needed[0],
        num_bins[1] + bin_needed[1],
        num_bins[2] + bin_needed[2],
        per_bin);
  }

  auto space_spec() const
  {
    // block the bins in each direction (will keep the fourth dimension on the
    // same unit)
    return dash::Pattern<4>(
        num_bins[0],
        num_bins[1],
        num_bins[2],
        per_bin + bin_buffer,
        dash::BLOCKED,
        dash::BLOCKED,
        dash::BLOCKED,
        dash::BLOCKED);
  }

  // number of atoms per bin (~<=8)
  auto per_bin_spec() const
  {
    // note this distribution means that the lookup is local for each unit's bins
    return dash::Pattern<3>(
        num_bins[0],
        num_bins[1],
        num_bins[2],
        dash::BLOCKED,
        dash::BLOCKED,
        dash::BLOCKED);
  }

  auto bin_spec() const
  {
    return space_spec();
  }

  uint32_t num_steps = 100;
  uint32_t size      = 32;
  uint32_t num_bins[3];
  float    bin_inv[3];
  uint32_t bin_needed[3];
  float    bin_size[3];
  int      per_bin    = 8;
  int      bin_buffer = 6;

  // input variables
  Input input;

  // simulation space specification
  std::valarray<md_type> box, boxlo, boxhi;
  md_type                volume    = 1.0;

  uint32_t               num_atoms = 1;
  md_type                mass      = 1.0;

  constexpr static md_type small  = 1.0e-6;
  constexpr static md_type factor = 0.999;

  friend std::ostream& operator<<(std::ostream& os, const Config& cfg)
  {
    os << "Config" << std::endl << std::endl;
    os << "\tnum_steps\t" << cfg.num_steps << std::endl;
    os << "\tsize\t\t" << cfg.size << std::endl;
    os << "\tnum_bins\t(" << cfg.num_bins[0] << ", " << cfg.num_bins[1]
       << ", " << cfg.num_bins[2] << ")" << std::endl;
    os << "\tbin_inv\t\t(" << cfg.bin_inv[0] << ", " << cfg.bin_inv[1] << ", "
       << cfg.bin_inv[2] << ")" << std::endl;
    os << "\tbin_needed\t(" << cfg.bin_needed[0] << ", " << cfg.bin_needed[1]
       << ", " << cfg.bin_needed[2] << ")" << std::endl;
    os << "\tbin_size\t(" << cfg.bin_size[0] << ", " << cfg.bin_size[1]
       << ", " << cfg.bin_size[2] << ")" << std::endl;
    os << "\tbox\t\t(" << cfg.box[0] << ", " << cfg.box[1] << ", "
       << cfg.box[2] << ")" << std::endl;
    os << "\tboxlo\t\t(" << cfg.boxlo[0] << ", " << cfg.boxlo[1] << ", "
       << cfg.boxlo[2] << ")" << std::endl;
    os << "\tboxhi\t\t(" << cfg.boxhi[0] << ", " << cfg.boxhi[1] << ", "
       << cfg.boxhi[2] << ")" << std::endl;
    os << "\tvolume\t\t" << cfg.volume << std::endl;
    os << "\tnum_atoms\t" << cfg.num_atoms << std::endl;
    os << "\tmass\t\t" << cfg.mass << std::endl;
    return os;
  }
};

// RNG
const uint32_t IA   = 16807;
const uint32_t IM   = 2147483647;
const double   AM   = 1.0 / IM;
const uint32_t IQ   = 127773;
const uint32_t IR   = 2836;
const uint32_t MASK = 123459876;
