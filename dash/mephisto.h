#pragma once

#include <mephisto/entity>
#include <mephisto/execution>
#include <libdash.h>
#include <patterns/local_pattern.h>


auto const Dim = 1;

class MephistoConfig {
  using SizeT   = dash::default_index_t;
  using EntityT = mephisto::Entity<Dim, SizeT, alpaka::acc::AccCpuSerial>;
  using Queue   = alpaka::queue::QueueCpuSync;
  using Context = mephisto::execution::AlpakaExecutionContext<EntityT, Queue>;
  using BasePattern = dash::BlockPattern<Dim>;
  using PatternT    = patterns::BalancedLocalPattern<BasePattern, EntityT>;
  using ExecutorT   = mephisto::execution::AlpakaExecutor<Context>;

private:
  Context _context;
  ExecutorT _executor;

public:
  MephistoConfig() : _context(), _executor(&_context) {}

  Context &context() {
    return _context;
  }

  auto executor()
  {
    return _executor;
  }

  auto policy()
  {
    return mephisto::execution::make_parallel_policy(_executor);
  }
};
