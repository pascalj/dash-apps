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

  Atom(Velocity &v) : v(v) {}

  explicit Atom(std::valarray<double> nv)
    : v(nv)
  {
  }
  explicit Atom(std::valarray<double> nv, std::valarray<double> nf)
    : v(nv)
    , f(nf)
  {
  }

  Atom() = default;

  void unset() {
    v = Velocity();
    f = Float3D();
  }

  void pbc(double box_x, double box_y, double box_z) {
    v.x += v.x < 0 ? box_x : 0;
    v.y += v.y < 0 ? box_y : 0;
    v.z += v.z < 0 ? box_z : 0;
    v.x -= v.x >= box_x ? box_x : 0;
    v.y -= v.y >= box_y ? box_y : 0;
    v.z -= v.z >= box_z ? box_z : 0;
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
    temperatures.allocate(dash::size());
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

      if ((curCoord.sum() % config.per_bin == 2) && within_bounds) {
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
          auto bin = coord2bin(temp);
      printf(
          "%d, %d, %d, -> %f, %f, %f -> %zu %zu %zu\n",
          curCoord[0],
          curCoord[1],
          curCoord[2],
          temp[0],
          temp[1],
          temp[2],
          bin[0],
          bin[1],
          bin[2]
          );
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


  void pbc() {
    dash::for_each(
        atoms.begin(),
        atoms.end(),
        [x = config.box[0], y = config.box[1], z = config.box[2]](Atom& a) {
          a.pbc(x, y, z);
        });
  }

  void add_atom(
      const Atom&                                a,
      const std::valarray<double>                v,
      const std::valarray<dash::default_index_t> bin)
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
    auto lbegin =  atoms.lbegin();
    const auto lend =  atoms.lend();
    Velocity accu;

    for(;lbegin < lend; lbegin++) {
      const auto vel = lbegin->v;
      accu.x += vel.x;
      accu.y += vel.y;
      accu.z += vel.z;
    }

    /* auto const total_velocity = dash::reduce( */
    /*     atoms.begin(), */
    /*     atoms.end(), */
    /*     Velocity(), */
    /*     [](const Velocity& lv, const Atom& rhs) { */
    /*       const auto rv = rhs.v; */
    /*       return Velocity{lv.x + rv.x, lv.y + rv.y, lv.z + rv.z}; */
    /*     }); */

    dash::Array<Velocity> local_sums{per_bin.team().size()};
    local_sums[dash::myid()] = accu;

    local_sums.barrier();

    Velocity total_velocity = std::accumulate(
        local_sums.begin(),
        local_sums.end(),
        Velocity(),
        [](const Velocity& accu, auto globref) {
          const Velocity v = globref;
          return Velocity{accu.x + v.x, accu.y + v.y, accu.z + v.z};
        });
    std::cout << "total velocity: " << total_velocity << std::endl;

    Velocity avg{total_velocity.x / config.num_atoms,
                 total_velocity.y / config.num_atoms,
                 total_velocity.z / config.num_atoms};


    dash::for_each(atoms.begin(), atoms.end(), [avg] (Atom &a) {
        a.v.x = -avg.x;
        a.v.y = -avg.y;
        a.v.z = -avg.z;
    });

    dash::barrier();
    const auto factor = std::sqrt(config.input.initial_temp / temperature());

    std::cout << "temp factor" << factor << std::endl;

    dash::for_each(atoms.begin(), atoms.end(), [avg, factor](Atom& a) {
      a.v.x *= factor;
      a.v.y *= factor;
      a.v.z *= factor;
    });

    /* dash::for_each_with_index( */
    /*     per_bin.begin(), per_bin.end(), [=](int bin_count, size_t bin_index) { */
    /*       auto const coords = per_bin.pattern().coords(bin_index); */
    /*       for (uint8_t i = 0; i < bin_count; i++) { */
    /*         atoms[coords[0]][coords[1]][coords[2]][i]; */
    /*       } */
    /*     }); */

  }

  double temperature()
  {
    double temp = 0;
    double local_temp = std::accumulate(
        atoms.lbegin(),
        atoms.lend(),
        0.0,
        [mass = config.mass](double accu, const Atom& a) {
          const auto v = a.v;
          return (v.x * v.x + v.y * v.y + v.z * v.z) * mass;
        });
    temperatures[dash::myid()] = local_temp;
    double global_temp =
        dash::reduce(temperatures.begin(), temperatures.end(), 0.0);
    return global_temp * config.input.temp_scale;
  }

  void bin_atoms()
  {
    // for each bin
    dash::for_each_with_index(
        per_bin.begin(),
        per_bin.end(),
        [=](int& bin_count, size_t bin_index) {
          auto const coords = per_bin.pattern().coords(bin_index);
          // for each atom in the bin
          for (uint8_t i = 0; i < bin_count; i++) {
            auto           a     = atoms[coords[0]][coords[1]][coords[2]][i];
            const Float3D& a_pos = pos[coords[0]][coords[1]][coords[2]][i];
            const auto     dest_bin = coord2bin({a_pos.x, a_pos.y, a_pos.z});
            auto const     val_coords =
                std::valarray<dash::default_index_t>(coords.data(), 3);
            if ((dest_bin != val_coords).min()) {
              // move the atom iff its bin changed
              add_atom(a, {a_pos.x, a_pos.y, a_pos.z}, dest_bin);

              if (i != bin_count - 1) {
                atoms[coords[0]][coords[1]][coords[2]][i] = a;
              }
              a = Atom();
              // TODO: atomic
              bin_count--;
            }
          }
        });
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

  std::valarray<dash::default_index_t> coord2bin(const std::valarray<md_type> x) const
  {
    std::valarray<dash::default_index_t> temp(
        static_cast<dash::default_index_t>(0), 3);

    for(uint8_t i = 0; i < 3; i++) {
      const float mask = x[i] >= config.box[i] ? 1 : 0;
      temp[i] = ((x[i] - config.box[i] * mask) * config.bin_inv[i]);
    }

    return temp;
  }

  const Config& config;
  dash::NArray<Float3D, 4> pos;
  dash::NArray<Float3D, 4> space;
  // actual atoms (x*y*z*per_bin[x][y][z])
  dash::NArray<Atom, 4> atoms;
  // Real positions of atoms
  // number of atoms per bin
  dash::NArray<int, 3> per_bin;
  dash::Array<double>  temperatures;
};


