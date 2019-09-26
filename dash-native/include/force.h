#pragma once

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

  void compute() {
    eng_vdwl = 0;
    virial = 0;
  }
};
