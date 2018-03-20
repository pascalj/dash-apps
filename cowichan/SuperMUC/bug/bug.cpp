#include <libdash.h>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::cin;

using value       = struct{ uint row, col; };

using pattern_t   = dash::CSRPattern< 1, dash::ROW_MAJOR, long >;
using extent_t    = pattern_t::size_type;
using res_array_t = dash::Array<value, long, pattern_t>;

static int myid;


int main( int argc, char* argv[] )
{
  dash::init( &argc,&argv );
  myid = dash::myid( );

  
  std::vector<extent_t> local_sizes_source;
  
  int i;  
  
  for( i=0; i < dash::Team::All().size(); ++i ){
    local_sizes_source.push_back(20+i);
  }

  pattern_t pattern_source( local_sizes_source );
  dash::Array<value, long, pattern_t> source ( pattern_source );

  i=0;  
  
  for( value * src = source.lbegin(); src < source.lend(); ++src, ++i )
  {
    src->col = i * myid;
    src->row = i + myid;
  }
  
  if(0==myid) cout << "size of source:" << source.size() << endl;
  value * target = new value[ source.size() ];

  if(0==myid) cout << "before dash::copy" << endl;
  
  dash::copy(source.begin(), source.end(), target);
  
  if(0==myid) cout << "after dash::copy" << endl;
  
  
  dash::finalize( );
}