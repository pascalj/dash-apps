#pragma once

#include <atom.h>
#include <init.h>
#include <algorithm>

struct Neighbor {
  // reference to an atom, either local or in a mirrored bin
  Atom*               a;
  dash::GlobPtr<Atom, dash::GlobStaticMem<dash::HostSpace>> glob_ptr;

  void update() {
    if (glob_ptr != nullptr) {
      *a = *glob_ptr;
    }
  }
};

struct Neighbors {
  Neighbors() = delete;

  Neighbors(const Config& cfg)
    : config(cfg)
  {
  }

  void rebuild(Atoms& atoms)
  {
    // reserve local neighbor lists
    neighs.clear();
    neighs.resize(atoms.atoms.local_size());
    const auto  nsize   = neighs.size();
    const auto  myid    = atoms.per_bin.team().myid();
    const auto& bin_pat = atoms.per_bin.pattern();

    const auto local_blocks = bin_pat.local_blockspec().size();
    // iterate over all local blocks
    const auto lindex_range =
        dash::local_index_range(atoms.atoms.begin(), atoms.atoms.end());
    const auto lbase = (atoms.atoms.begin() + lindex_range.begin).local();

    for (int b = 0; b < local_blocks; b++) {
      const auto offsets = bin_pat.local_block(b).offsets();
      const auto extents = bin_pat.local_block(b).extents();
      const auto xstart  = offsets[0];
      const auto ystart  = offsets[1];
      const auto zstart  = offsets[2];
      const auto xend    = extents[0];
      const auto yend    = extents[1];
      const auto zend    = extents[2];

      for (uint32_t x = xstart; x < xend; x++) {
        for (uint32_t y = ystart; y < yend; y++) {
          for (uint32_t z = zstart; z < zend; z++) {
            const int  bin_size = atoms.per_bin[x][y][z];
            const auto l_bin_offset =
                atoms.atoms.pattern().local_at({x, y, z, 0});
            const auto l_bin = lbase + l_bin_offset;

            neighbor_bins(
                atoms,
                x,
                y,
                z,
                [&, this](Atom* n_bin, auto glob_ptr, size_t nbin_size) {
                  // for each atom in the bin
                  for (uint32_t j{0}; j < bin_size; j++) {
                    // for each neighbor in the neighbor bin
                    for (uint32_t i{0}; i < nbin_size; i++) {
                      if (i == j && n_bin == l_bin) {
                        continue;
                      }
                      auto a_pos = l_bin[j].pos;
                      auto n_pos = n_bin[i].pos;

                      auto rsq = distance(a_pos, n_pos);
                      if (rsq <= config.input.cutneighsq) {
                        auto update_ptr =
                            glob_ptr == nullptr ? glob_ptr : glob_ptr + i;
                        Neighbor n{n_bin + i, update_ptr};
                        neighs[l_bin_offset + j].push_back(n);
                      }
                    }
                  }
                });
          }
        }
      }
    }
  }

  void update_positions(Atoms& atoms)
  {
    for (auto& atom_neighs : neighs) {
      for (auto& neigh : atom_neighs) {
        neigh.update();
      }
    }
  }

  template <typename F>
  void neighbor_bins(Atoms& a, uint32_t x, uint32_t y, uint32_t z, F&& f)
  {
    using gptr_t = dash::GlobPtr<Atom, dash::GlobStaticMem<dash::HostSpace>>;
    const auto minx = std::max(uint32_t{0}, x - config.bin_needed[0]);
    const auto miny = std::max(uint32_t{0}, y - config.bin_needed[1]);
    const auto minz = std::max(uint32_t{0}, z - config.bin_needed[2]);
    const auto maxx = std::min(config.num_bins[0], x + config.bin_needed[0]);
    const auto maxy = std::min(config.num_bins[1], y + config.bin_needed[1]);
    const auto maxz = std::min(config.num_bins[2], z + config.bin_needed[2]);

    const auto& pat = a.atoms.pattern();
    for (auto nx = minx; nx < maxx; nx++) {
      for (auto ny = miny; ny < maxy; ny++) {
        for (auto nz = minz; nz < maxz; nz++) {
          const auto gindex    = pat.global_at({nx, ny, nz, 0});
          const int  nbin_size = a.per_bin[nx][ny][nz];

          Atom*  bin_ptr;
          gptr_t ptr;
          if (pat.is_local(gindex)) {
            bin_ptr = (a.atoms.begin() + gindex).local();
          }
          else {
            // grow mirrored atoms list
            const auto insert_pos = mirrored_atoms.size();
            mirrored_atoms.resize(insert_pos + nbin_size);
            bin_ptr = mirrored_atoms.data() + insert_pos;
            ptr = gptr_t{a.atoms.begin() + gindex};
            std::cout << "ptr: " << a.atoms.begin() << std::endl;
            dash::copy(ptr, ptr + nbin_size, bin_ptr);
            std::cout << "last: " << mirrored_atoms.back() << std::endl;
            std::cout << "size: " << mirrored_atoms.size() << std::endl;
          }
          f(bin_ptr, ptr, nbin_size);
        }
      }
    }
  }

  /* template<typename CoordsT> */
  /* auto& neighbors_for(CoordsT coords) { */
  /*   return neighs[gindex]; */
  /* } */

  std::vector<std::vector<Neighbor>> neighs;
  std::vector<Atom>                  mirrored_atoms;
  const Config&                      config;
};

