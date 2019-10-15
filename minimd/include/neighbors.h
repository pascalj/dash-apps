#pragma once

#include <atom.h>
#include <init.h>
#include <algorithm>
#include <map>

struct Neighbor {
  using ref_t = dash::GlobAsyncRef<Atom>;
  // reference to an atom, either local or in a mirrored bin
  Atom*               a;
  ref_t glob_ptr;

  Neighbor(Atom* a, ref_t&& ref)
    : a(a)
    , glob_ptr(std::move(ref))
  {
  }

  void update() {
    if (glob_ptr.dart_gptr() != DART_GPTR_NULL) {
      glob_ptr.get(a);
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
    const auto& bin_pat = atoms.per_bin.pattern();

    const auto local_blocks = bin_pat.local_blockspec().size();
    // iterate over all local blocks
    const auto lindex_range =
        dash::local_index_range(atoms.atoms.begin(), atoms.atoms.end());
    const auto lbase = (atoms.atoms.begin() + lindex_range.begin).local();

    for (size_t b = 0; b < local_blocks; b++) {
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
            const uint32_t bin_size = atoms.per_bin[x][y][z];
            const auto     l_bin_offset =
                atoms.atoms.pattern().local_at({x, y, z, 0});
            const auto l_bin = lbase + l_bin_offset;

            neighbor_bins(
                atoms,
                x,
                y,
                z,
                [&, this](Atom* n_bin, auto glob_ptr, uint32_t nbin_size) {
                  // for each atom in the bin
                  for (uint32_t j = 0; j < bin_size; j++) {
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
                        dash::GlobAsyncRef<Atom> ref(update_ptr.dart_gptr());
                        neighs[l_bin_offset + j].emplace_back(n_bin + i, std::move(ref));
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
    // Flush all updates
    auto gptr = atoms.per_bin.begin().dart_gptr();
    dart_flush_all(gptr);
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
          const auto gindex     = pat.global_at({nx, ny, nz, 0});
          const auto nbin_index = a.per_bin.pattern().global_at({nx, ny, nz});
          const int  nbin_size  = a.per_bin[nx][ny][nz];

          Atom*  bin_ptr;
          gptr_t ptr;
          if (pat.is_local(gindex)) {
            assert(a.per_bin.pattern().is_local(nbin_index));
            bin_ptr = (a.atoms.begin() + gindex).local();
          }
          else {
            assert(!a.per_bin.pattern().is_local(nbin_index));
            auto iter = a.atoms.begin() + gindex;
            assert(!iter.is_local());
            ptr       = static_cast<gptr_t>(iter);
            // grow mirrored atoms list
            if (!mirrored_bins.count(nbin_index)) {
              auto mirrored_bin =
                  mirrored_bins.emplace(nbin_index, nbin_size);
              dash::copy(
                  iter, iter + nbin_size, mirrored_bin.first->second.data());
            }
            bin_ptr = mirrored_bins[nbin_index].data();
          }
          f(bin_ptr, ptr, nbin_size);
        }
      }
    }
  }

  std::vector<std::vector<Neighbor>>   neighs;
  // simple store for mirrored bins
  std::map<index_t, std::vector<Atom>> mirrored_bins;
  const Config&                        config;
};

