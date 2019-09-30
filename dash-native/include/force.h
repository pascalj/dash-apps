#pragma once

#include "atom.h"

struct Force {
  float cutforcesq;
  float eng_vdwl;

  float virial;
  float mass;
};

struct ForceLJ : public Force {
  ForceLJ() = delete;

  ForceLJ(const float cf) : cutforcesq(cf * cf) {}

  float cutforcesq;

  void compute(Atoms& atoms) {
    eng_vdwl = 0;
    virial = 0;

    auto& per_bin = atoms.per_bin;
    auto& bins = atoms.atoms;

    dash::for_each(per_bin.begin(), per_bin.end(), [](Atom& a) {
        a.v = Velocity();
        a.f = Float3D();
    });

    dash::for_each(per_bin.begin(), per_bin.end(), [](Atom& a) {
    });
  }
};
