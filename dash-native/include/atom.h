#pragma once

#include <libdash.h>
#include <algorithm>
#include <valarray>
#include <vector>
#include <tuple>
#include "init.h"

using Float3D = std::array<double, 3>;
using Velocity = Float3D;
using v3int = std::valarray<int>;
using v3 = std::valarray<double>;
using pos_t = std::array<dash::default_index_t, 4>;

struct Atom {
  Velocity v;
  Float3D  f{{0.0, 0.0, 0.0}};
  Float3D  pos;

  explicit Atom(Float3D nv)
    : v(nv)
  {
  }
  explicit Atom(Float3D nv, Float3D nf)
    : v(nv)
    , f(nf)
  {
  }
  explicit Atom(Float3D nv, Float3D nf, Float3D npos)
    : v(nv)
    , f(nf)
    , pos(npos)
  {
  }

  Atom() = default;

  void unset() {
    v = Velocity();
    f = Float3D();
  }


  /**
   * Reset boundary conditions
   */
  void pbc(double box_x, double box_y, double box_z) {
    v[0] += v[0] < 0 ? box_x : 0;
    v[1] += v[1] < 0 ? box_y : 0;
    v[2] += v[2] < 0 ? box_z : 0;
    v[0] -= v[0] >= box_x ? box_x : 0;
    v[1] -= v[1] >= box_y ? box_y : 0;
    v[2] -= v[2] >= box_z ? box_z : 0;
  }

  friend std::ostream& operator<<(std::ostream& os, const Atom& a)
  {
    os << "Atom([" << a.v[0] << "," << a.v[1] << "," << a.v[2] << "], ";
    os << "[" << a.f[0] << "," << a.f[1] << "," << a.f[2] << "], ";
    os << "[" << a.pos[0] << "," << a.pos[1] << "," << a.pos[2] << "])";
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
    space.allocate(cfg.space_spec());
    // bins x*y*z*8
    atoms.allocate(cfg.bin_spec());
    per_bin.allocate(cfg.per_bin_spec());
    temperatures.allocate(dash::size());
    t_eng_arr.allocate(dash::size());
    t_vir_arr.allocate(dash::size());
  }

  template <typename IndexT>
  auto get_atom(const std::array<IndexT, 3>& coords, const int offset_in_bin)
  {
    return atoms[coords[0]][coords[1]][coords[2]][offset_in_bin];
  }

  template <typename IndexT>
  auto get_atom(const std::array<IndexT, 4>& coords)
  {
    return atoms[coords[0]][coords[1]][coords[2]][coords[3]];
  }

  template <typename IndexT>
  auto get_pos(const std::array<IndexT, 3>& coords, const IndexT i)
  {
    return atoms[coords[0]][coords[1]][coords[2]][i].pos;
  }

  template <typename IndexT>
  auto get_pos(const std::array<IndexT, 4>& coords)
  {
    return atoms(coords).get().pos;
  }

  /* template <typename IndexT> */
  /* pos_t get_neighbor(const std::array<IndexT, 3>& coords, const int offset) */
  /* { */
  /*   return neighbors[coords[0]][coords[1]][coords[2]][offset]; */
  /* } */

  template <typename F>
  void each_atom_with_offset(F&& f)
  {
    dash::for_each_with_index(
        per_bin.begin(),
        per_bin.end(),
        [f, p = per_bin.pattern(), this](int natoms, size_t bin_index) {
          auto const coords = p.coords(bin_index);

          std::array<dash::default_index_t, 4> a_begin = {
              coords[0], coords[1], coords[2], 0};
          auto begin =
              dash::local(atoms.begin() + atoms.pattern().global_at(a_begin));
          for (int i = 0; i < natoms; i++) {
            a_begin[3]++; 
            f(*(begin + i), a_begin);
          }
        });
  }

  void initial_integrate()
  {
    each_atom_with_offset(
        [this](
            Atom&                                       a,
            const std::array<dash::default_index_t, 4>& coords) {
          for (int i = 0; i < 3; i++) {
            a.v[i] += a.f[i] * config.input.dtforce;
            a.pos[i] += a.v[i] * config.input.dt;
          }
        });
  }

  void final_integrate() {
    dash::for_each(atoms.begin(), atoms.end(), [this] (Atom& a) {
        for(int i = 0; i < 3; i++) {
          a.v[i] += config.input.dtforce * a.f[i];
        }
    });
  }


  /**
   * Create all needed atoms from a PRNG distributed across box, stored into
   * bins.
   */
  void create_atoms()
  {
    std::valarray<int> lo{0, 3}, hi{0, 3};

    for (uint8_t i = 0; i < 3; i++) {
      lo[i] = std::max(0.0, config.boxlo[i] / (0.5 * config.input.lattice));
      hi[i] = std::min(
          2.0 * config.input.problem_size[i] - 1,
          config.boxhi[i] / (0.5 * config.input.lattice));
    }

    std::valarray<int>    o(0, 3), s(0, 3), curCoord(0, 3);
    std::valarray<double> temp(0.0, 3);
    Float3D v;
    size_t                total = 0;

    while (o[2] * config.per_bin <= hi[2]) {
      curCoord = o * config.per_bin + s;

      bool within_bounds = (curCoord >= lo).min() && (curCoord <= hi).min();

      if ((curCoord.sum() % config.per_bin == 2) && within_bounds) {
        for (uint8_t i = 0; i < 3; i++) {
          temp[i] = curCoord[i] * 0.5 * config.input.lattice;
        }

        within_bounds =
            (temp >= config.boxlo).min() && (temp < config.boxhi).min();

        const Input& in = config.input;
        if (within_bounds) {
          int n = curCoord[2] * (2 * in.problem_size[1]) *
                      (2 * in.problem_size[0]) +
                  curCoord[1] * (2 * in.problem_size[0]) + curCoord[0] + 1;
          for (uint8_t i = 0; i < 3; i++) {
            for (uint8_t m = 0; m < 5; m++) {
              pmrand(n);
            }
            v[i] = pmrand(n);
          }
          auto bin = coord2bin(temp);
          Atom a(v);
          for (uint8_t i = 0; i < 3; i++) {
            a.pos[i] = temp[i];
          }
          add_atom(a, temp, coord2bin(temp));
          total++;
        }
      }

      s[0]++;

      if (s[0] == config.per_bin) {
        s[0] = 0;
        s[1]++;
      }
      if (s[1] == config.per_bin) {
        s[1] = 0;
        s[2]++;
      }
      if (s[2] == config.per_bin) {
        s[2] = 0;
        o[0]++;
      }
      if (o[0] * config.per_bin > hi[0]) {
        o[0] = 0;
        o[1]++;
      }
      if (o[1] * config.per_bin > hi[1]) {
        o[1] = 0;
        o[2]++;
      }
    }

    /* std::vector<int> bin_hist(config.per_bin + config.bin_buffer); */
    /* int              i = 0; */
    /* for (int count : per_bin) { */
    /*   bin_hist[count]++; */
    /*   std::cout << " " << std::setw(2) << count; */
    /*   if (i++ % 26 == 0) { */
    /*     std::cout << std::endl; */
    /*   } */
    /* } */
    /* std::cout << "hist: "; */
    /* for (auto bin : bin_hist) { */
    /*   std::cout << " " << bin; */
    /* } */
    /* std::cout << std::endl; */
    /* ; */
  }

  void pbc()
  {
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
    int bin_atoms = per_bin[bin[0]][bin[1]][bin[2]]++;
    if (bin_atoms >= config.per_bin + config.bin_buffer) {
      std::cout << "max atoms exhausted" << bin_atoms << std::endl;
    }
    atoms[bin[0]][bin[1]][bin[2]][bin_atoms] = a;
  }

  void create_velocity() {
    auto       lbegin = atoms.lbegin();
    const auto lend   = atoms.lend();
    Velocity   accu;

    for(;lbegin < lend; lbegin++) {
      const auto vel = lbegin->v;
      accu[0] += vel[0];
      accu[1] += vel[1];
      accu[2] += vel[2];
    }

    /* auto const total_velocity = dash::reduce( */
    /*     atoms.begin(), */
    /*     atoms.end(), */
    /*     Velocity(), */
    /*     [](const Velocity& lv, const Atom& rhs) { */
    /*       const auto rv = rhs.v; */
    /*       return Velocity{lv[0] + rv[0], lv[1] + rv[1], lv[2] + rv[2]}; */
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
          return Velocity{accu[0] + v[0], accu[1] + v[1], accu[2] + v[2]};
        });

    Velocity avg{total_velocity[0] / config.num_atoms,
                 total_velocity[1] / config.num_atoms,
                 total_velocity[2] / config.num_atoms};

    dash::for_each_with_index(
        atoms.begin(),
        atoms.end(),
        [avg, this](Atom& a, dash::default_index_t index) {
          auto c = atoms.pattern().coords(index);
          if (per_bin[c[0]][c[1]][c[2]] > c[3]) {
            a.v[0] -= -avg[0];
            a.v[1] -= -avg[1];
            a.v[2] -= -avg[2];
          }
        });

    dash::barrier();
    const auto factor = std::sqrt(config.input.initial_temp / temperature());

    dash::for_each_with_index(
        atoms.begin(),
        atoms.end(),
        [avg, this, factor](Atom& a, index_t index) {
          auto c = atoms.pattern().coords(index);
          if (per_bin[c[0]][c[1]][c[2]] > c[3]) {
            a.v[0] *= factor;
            a.v[1] *= factor;
            a.v[2] *= factor;
          }
        });
  }

  double temperature()
  {
    double local_temp = std::accumulate(
        atoms.lbegin(),
        atoms.lend(),
        0.0,
        [mass = config.mass](double accu, const Atom& a) {
          const auto v = a.v;
          return accu + (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * mass;
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
        [=](uint32_t& bin_count, size_t bin_index) {
          auto const coords = per_bin.pattern().coords(bin_index);
          // for each atom in the bin
          for (uint8_t i = 0; i < bin_count; i++) {
            Atom           a     = atoms[coords[0]][coords[1]][coords[2]][i];
            const Float3D& a_pos = a.pos;
            const auto     dest_bin = coord2bin({a_pos[0], a_pos[1], a_pos[2]});
            auto const     val_coords =
                std::valarray<dash::default_index_t>(coords.data(), 3);
            if ((dest_bin != val_coords).min()) {
              // move the atom iff its bin changed
              add_atom(a, {a_pos[0], a_pos[1], a_pos[2]}, dest_bin);

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

  const Config&             config;
  dash::NArray<Float3D, 4>  space;
  dash::NArray<Atom, 4>     atoms;    // actual atoms (x*y*z*per_bin[x][y][z])
  dash::NArray<uint32_t, 3> per_bin;  // number of atoms per bin
  dash::Array<double>       temperatures;  // temperature on each unit
  dash::Array<double>       t_eng_arr;     // temperature on each unit
  dash::Array<double>       t_vir_arr;     // temperature on each unit
};

template <typename Arr>
auto distance(const Arr a, const Arr b)
{
  using value_t = typename Arr::value_type;
  value_t sum{0};
  for (size_t i = 0; i < std::tuple_size<Arr>::value; i++) {
    const auto ai = a[i];
    const auto bi = b[i];
    sum += (ai - bi) * (ai - bi);
  }
  return sum;
}

