#pragma once

#include "atom.h"
#include "neighbors.h"

template <typename T>
auto dot(const T& a, const T& b)
{
  constexpr size_t size = std::tuple_size<T>::value;
  typename T::value_type sum{0};
  for (size_t i = 0; i < size; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

struct Force {
  float eng_vdwl;

  float virial;
  float mass;
};

struct ForceLJ : public Force {
  ForceLJ() = delete;
  float cutforcesq;

  ForceLJ(const float cf) : cutforcesq(cf * cf) {}

  void compute(Atoms& atoms, Neighbors &neighs, bool store = false) {
    eng_vdwl = 0;
    virial = 0;

    auto& per_bin = atoms.per_bin;
    auto& bins = atoms.atoms;

    dash::for_each(atoms.atoms.begin(), atoms.atoms.end(), [](Atom& a) {
      a.f = Float3D();
    });

    const auto l_begin =
        dash::local_index_range(atoms.per_bin.begin(), atoms.per_bin.end()).begin;

    dash::for_each_with_index(
        per_bin.begin(),
        per_bin.end(),
        [&, this](int natoms, size_t bin_index) {
          auto lbin_index = bin_index - l_begin;
          auto coords     = per_bin.pattern().coords(bin_index);

          assert(lbin_index >= 0);
          std::array<md_type, 3> del;

          auto a_begin =
              atoms.atoms[coords[0]][coords[1]][coords[2]].begin().local();
          const auto atom_index = atoms.atoms.pattern().local_at(
              {coords[0], coords[1], coords[2], 0});

          for (dash::default_index_t i = 0; i < natoms; i++) {
            auto atom = a_begin + i;
            for (const auto& neigh : neighs.neighs[atom_index + i]) {
              for(int d = 0; d < 3; d++) {
                del[d] = atom->pos[d] - neigh.a->pos[d];
              }

              auto rsq = dot(del, del);
              if (rsq < cutforcesq) {
                const auto sr2   = 1 / rsq;
                const auto sr6   = std::pow(sr2, 3);
                const auto force = 48.0 * sr6 * (sr6 - 0.5) * sr2;
                for (int b = 0; b < 3; b++) {
                  atom->f[b] += del[b] * force;
                }
              }
            }
          }
        });
  }
};
