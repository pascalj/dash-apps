#pragma once

#include <libdash.h>
#include <algorithm>
#include <valarray>
#include "init.h"

using Velocity = Float3D;
using v3int = std::valarray<int>;
using v3 = std::valarray<double>;

struct Atom {
  Velocity v;
  Float3D  f{{0.0, 0.0, 0.0}};

  Atom() = delete;

  explicit Atom(std::valarray<double> nv)
    : v(nv)
  {
  }
  explicit Atom(std::valarray<double> nv, std::valarray<double> nf)
    : v(nv)
    , f(nf)
  {
  }

  friend std::ostream& operator<<(std::ostream& os, const Atom& a)
  {
    os << "Atom([" << a.v.x << "," << a.v.y << "," << a.v.z << "], ";
    os << "[" << a.f.x << "," << a.f.y << "," << a.f.z << "])";
    return os;
  }
};

/**
 * The total collection of all atoms
 */
struct Atoms {
  Atoms() = delete;

  Atoms(const Config& cfg)
    : config(cfg)
  {
    pos.allocate(cfg.pos_spec());
    space.allocate(cfg.space_spec());
    // bins x*y*z*8
    atoms.allocate(cfg.bin_spec());
    per_bin.allocate(cfg.per_bin_spec());
  }

  void create_atoms() {
    std::valarray<int> lo{0, 3}, hi{0, 3};

    for(uint8_t i = 0; i < 3; i++) {
      lo[i] = std::max(0.0, config.boxlo[i] / (0.5 * config.input.lattice));
      hi[i] = std::min(
          2.0 * config.input.problem_size[i] - 1,
          config.boxhi[i] / (0.5 * config.input.lattice));
    }

    std::valarray<int>   o(0, 3), s(0, 3), curCoord(0, 3);
    std::valarray<double> temp(0.0, 3), v(0.0, 3);
    size_t total = 0;

    while(o[2] * config.per_bin <= hi[2]) {
      curCoord = o * config.per_bin + s;

      bool within_bounds = (curCoord >= lo).min() && (curCoord <= hi).min();

      if ((curCoord.sum() % config.per_bin == 0) && within_bounds) {
        for (uint8_t i = 0; i < 3; i++) {
          temp[i] = curCoord[i] * 0.5 * config.input.lattice;
        }

        within_bounds = (temp >= config.boxlo).min() && (temp < config.boxhi).min();

        const Input& in = config.input;
        if (within_bounds) {
          int n = curCoord[2] * (2 * in.problem_size[1]) *
                   (2 * in.problem_size[0]) + curCoord[1] *
                   (2 * in.problem_size[0]) + curCoord[0] + 1;
          for(uint8_t i = 0; i < 3; i++) {
            for(uint8_t m = 0; m < 5; m++) {
              pmrand(n);
            }
            v[i] = pmrand(n);
          }
          Atom a{v};
          add_atom(a, temp, coord2bin(temp));
          total++;
        }
      }

      s[0]++;

      if(s[0] == config.per_bin) {
        s[0] = 0;
        s[1]++;
      }
      if(s[1] == config.per_bin) {
        s[1] = 0;
        s[2]++;
      }
      if(s[2] == config.per_bin) {
        s[2] = 0;
        o[0]++;
      }
      if(o[0] * config.per_bin > hi[0]) {
        o[0] = 0;
        o[1]++;
      }
      if(o[1] * config.per_bin > hi[1]) {
        o[1] = 0;
        o[2]++;
      }
    }
    std::cout << "total " << total << std::endl;

    std::vector<int> bin_hist(config.per_bin + config.bin_buffer);
    int i = 0;
    for(int count : per_bin) {
      bin_hist[count]++;
      std::cout << " " << std::setw(2) << count;
      if (i++ % 26 == 0) {
        std::cout << std::endl;
      }
    }
    std::cout << "hist: ";
    for(auto bin : bin_hist) {
      std::cout << " " << bin;
    }
    std::cout << std::endl;;
  }

  void add_atom(
      Atom& a, const std::valarray<double> v, const std::valarray<int> bin)
  {
    
    /* std::cout << "adding atom: " << a << " -> " */
    /*           << "(" << bin[0] << "," << bin[1] << "," << bin[2] << ")" */
    /*           << std::endl; */
    int bin_atoms = per_bin[bin[0]][bin[1]][bin[2]]++;
    if(bin_atoms >= config.per_bin + config.bin_buffer) {
      std::cout << "max atoms exhausted" << bin_atoms << std::endl;
    }
    atoms[bin[0]][bin[1]][bin[2]][bin_atoms] = a;

  }

  void create_velocity() {
    
  }

  double pmrand(int& n) const
  {
    int    k;
    double ans;

    k = n / IQ;
    n = IA * (n - k * IQ) - IR * k;
    if (n < 0) {
      n += IM;
    }
    ans = AM * n;
    return ans;
  }

  std::valarray<int> coord2bin(const std::valarray<md_type> x) const
  {
    v3int cur(0, 3), temp(0, 3);

    for(uint8_t i = 0; i < 3; i++) {
      const float mask = x[i] >= config.box[i] ? 1 : 0;
      temp[i] = ((x[i] - config.box[i] * mask) * config.bin_inv[i]);
      /* std::cout << " (" << x[i] << ", " << temp[i] << ")"; */
    }

    return temp;
  }

  const Config& config;
  dash::NArray<Float3D, 4> pos;
  dash::NArray<Float3D, 4> space;
  // actual atoms (x*y*z*per_bin[x][y][z])
  dash::NArray<Atom, 4> atoms;
  // number of atoms per bin
  dash::NArray<int, 3>     per_bin;
};


