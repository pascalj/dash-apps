// #ifndef DASH_ENABLE_LOGGING
// #define DASH_ENABLE_LOGGING
// #endif

#include <libdash.h>

#include "../Terminal_Color.h"
#include <cstdlib>

#define DEBUG

#ifdef DEBUG
  #include <chrono>
  #include <thread>
#endif

#define MAX_KEY               99
#define MIN_NUM_ELEM_PER_UNIT 10


using std::cout   ;
using std::cin    ;
using std::endl   ;
using std::vector ;
using std::pair   ;

using dash::Team   ;
using dash::Array  ;
using dash::NArray ;
using dash::Shared ;
using dash::fill   ;

using uint       = unsigned  int ;
using uchar      = unsigned char ;
using POI_T      =          int  ;  //this type musst be signed!
using MATRIX_T   =         uchar ;

using pattern_t  = dash::CSRPattern< 1, dash::ROW_MAJOR, int >;
using extent_t   = pattern_t::size_type;

using Point      = struct{ MATRIX_T value; uint row, col;   };
using pointRange = struct{ Point * curr, * end;             };  // unused until now
using unitRange  = struct{ MATRIX_T begin, end; uint count; };  // unitValueRange

// static variables
static struct InputPar { uint nrows, ncols; } in;
static uint   nelts;
static int    myid ;

template< typename T >
void print2d( const T& mat ) {
  for( int i = 0; i < mat.extent(0); i++ ) {
    for( int j = 0; j < mat.extent(1); j++ ) {
      cout << std::setw(3) << static_cast<uint>( mat(i,j) )<< " ";
    }
    cout << endl;
  }
}

bool operator<(const Point& lhs, const Point& rhs)
{
  return lhs.value < rhs.value;
}

std::ostream& operator<<(std::ostream& os, const Point& p)
{
  #ifdef DEBUG
    return os << "(" << fmt( p.value, FCYN ) << "," << fmt( p.row, FGREEN ) << "," << fmt( p.col, FGREEN ) << ")-";
  #else
    return os << p.row << " " << p.col << "\n";
  #endif
}

#ifdef DEBUG
  using std::this_thread::sleep_for;
  inline void __sleep( uint const baseDur = 0, uint const mult = 10 )
  {
     uint SLEEP_TIME__ = (myid + 1) * mult + baseDur;
      sleep_for(std::chrono::milliseconds(SLEEP_TIME__)); 
  }
#endif

/*
 * One unit has the job to read in the parameters.
 * Because there's always a unit0, it reads the input parameter and
 * distributes them to the rest of the units.
 */
inline void ReadRowsNCols( )
{
  Shared<InputPar> input_transfer;

  if(0 == myid)
  {
    cin >> in.nrows;
    cin >> in.ncols;

    input_transfer.set(in);
  }
  input_transfer.barrier();
  in = input_transfer.get();
}


template< typename T = MATRIX_T >
inline void ReadMatricesAndNelts( NArray<T,2>& randMat, NArray<bool,2>& threshMask )
{
  Shared<uint> nelts_transfer;

  if(0 == myid)
  {
    //read matrices
    T tmp;
      for ( auto i : randMat ){
        scanf( "%u", &tmp )  , i = tmp;
      }
      bool tmpB;
      for ( auto i : threshMask ){
        scanf( "%u", &tmpB ) , i = tmpB;
      }
    cin >> nelts;

    nelts_transfer.set(nelts);
  }
  nelts_transfer.barrier();
  nelts = nelts_transfer.get();
}


template<typename T = MATRIX_T, typename X = pair<POI_T, POI_T> >
inline void winnow(
              uint const   nrows      ,
              uint const   ncols      ,
    NArray<   T,2> const & randMat    ,
    NArray<bool,2> const & threshMask ,
              uint const   nelts      ,
         vector<X>       & result     )
{  
  Team & team   = dash::Team::All ( );
  size_t nUnits = team.size       ( );
  
  #ifdef DEBUG
    if( 0 == myid ) cout << "nUnits: " << nUnits << endl;
  #endif
  /* create global histo array for sorting
   * size += 1 for direct Index access -> histo[2]++ counts for value 2
   * size += 1 for additional value at the end of the 
   * histogram -> used for counter how many values were found
   */
  Array<uint> histo( (MAX_KEY + 2) * nUnits, dash::BLOCKED );
  // fill( histo.begin( ), histo.end( ), 0 );
  for(uint * it = histo.lbegin( ); it < histo.lend( ); ++it){*it = 0;}
  
  // local found points are gathered in this vector
  vector<Point> pointsLocal;
    
  // returns a object with the global row and column of the the local coordinates {0,0}
  auto globIndex = randMat.pattern( ).global( {0,0} );
  
  uint         gRow   = globIndex[0]           ;
  uint         gCol   = globIndex[1]           ;
  uint       & found  = histo.local[MAX_KEY+1] ;
     T const * matrEl = randMat.lbegin( )      ;
  
  
 /* read in local part of mask - matrix combination
  * and while doing that:
  * generate histogram and count the values (-> found)
  */
  for ( bool const * maskEl = threshMask.lbegin( );  maskEl < threshMask.lend( );  ++maskEl , ++matrEl, ++gCol )
  {
    if( gCol == ncols ) gCol = 0, ++gRow; // end of row -> next row in matrix/mask
    if( *maskEl       )
    {
      pointsLocal.push_back( Point{ *matrEl, gRow, gCol } );
      ++histo.local[*matrEl];
    }
  }
  found = pointsLocal.size();
  
  #ifdef DEBUG  // print points found local
    dash::barrier();  // only needed for better IO Output
    __sleep( );
    
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( auto it : pointsLocal){ cout <<   it; }
    cout << endl;
    
    __sleep(20);
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": found via pointsLocal:" << fmt( pointsLocal.size()   , FRED, 2 )  << "\n";
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": found via histogram  :" << fmt( found                , FRED, 2 )  << "\n";
    // cout << "#" << fmt( myid, FBLUE, 2 ) << ": histogram.lsize:"       << fmt(found - histo.lbegin(), FRED, 2 )  << endl;
  #endif

  // is here a barrier needed? -> to wait for completion on unit0?!
  if( 0 != myid ) {
    dash::transform<uint>(
      histo.lbegin     ( ) ,
      histo.lend       ( ) ,
      histo.begin      ( ) ,
      histo.begin      ( ) ,
      dash::plus<uint> ( ) );
  }
  
  
  // in this array will be the distribution info for creating buckets
  unitRange * distr = static_cast<unitRange*>(  std::malloc( nUnits * sizeof(unitRange) )  );
  unitRange * const distr_end = distr + nUnits;


  // unit 0 have to wait for rest to finish adding their values
  dash::barrier();  
  
  if( 0 == myid ) { // calculate bucket distribution

    /* calculate how much elements each unit should hold ideally
     * increment by one for safety garuantees in distribution
     * that is the assumption (ideal * nUnits > foundAllSize) == true
     */
     uint ideal = std::max( static_cast<size_t>( MIN_NUM_ELEM_PER_UNIT ), (found / nUnits) + 1 );
     
    #ifdef DEBUG
       __sleep();
      cout << "ideal number of elements per unit: " << fmt( ideal, FRED ) << endl;
      
      cout << "Histogram: ";
      for( size_t i = 0; i < histo.lsize(); ++i ) {
        if( histo.local[i] ) cout << fmt( i, FCYN ) << ":" << fmt( histo.local[i], FRED ) << ", ";
      } 
      cout << endl;
    #endif
     
  
    // begin - 1 for loop logic (starting with prefix ++)
    unitRange * uRPtr = distr - 1;
    
   /* the loop for calculation of the distribution got a bit more complex
    * because i wanted to iterate only once over "distr"
    * therefore no initialization beforehand is needed
    */ 
    if ( 1 == nUnits )
    {
      (++uRPtr)->begin = 0;
      uRPtr->end       = MAX_KEY;
      uRPtr->count     = found;
    
    }else{
      
      uint   acc    = 0;
         T   begin  = 0;
      uint * hisIt  = histo.lbegin();
      
      // actual calculation of distribution
      for( size_t i = 0; i < histo.lsize()-1 ; ++i ) {
        acc += *hisIt++;
        
        if( acc >= ideal ){
          
          (++uRPtr)->begin = begin;
          uRPtr->end = i;
          uRPtr->count = acc;
          
          begin = i+1;
          acc   = 0;
        }
      }
      uRPtr->count += acc;
      uRPtr->end = MAX_KEY;
    }
    // set the rest to zero
    while( ++uRPtr < distr_end ) { uRPtr->begin = 0; uRPtr->end = 0; uRPtr->count = 0; }
    
  } // end of unit 0 only part
  
  
  // convert own team unit ID into team unit ID with 0
  dash::team_unit_t TeamUnit0ID = team.myid( );
  TeamUnit0ID.id = 0;

  dart_ret_t ret = dart_bcast(
                      static_cast<void*>( distr )        ,  // buf 
                      nUnits * sizeof(unitRange)         ,  // nelem
                      DART_TYPE_BYTE                     ,  // dtype
                      TeamUnit0ID                        ,  // root
                      team.dart_id( )                    ); // team
  
  if( DART_OK != ret ) cout << "An Error while BCAST has occured!" << endl; 

  #ifdef DEBUG
    __sleep();
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( unitRange * i = distr; i < distr_end; ++i )
    {
        cout << "Range: " 
          << fmt( i->begin, FYEL ) 
          << "-" 
          << fmt( i->end, FGREEN )
          << " c:"
          << i->count
          << ", ";
      }
      cout << endl;
  #endif
  
 /* create dash::Array which will hold the data to be sorted.
  * therefore the CSRPattern of DASH will be used.
  * this allows different local sizes!
  */ 
  vector<extent_t> local_sizes;
  for( unitRange * uRPtr = distr; uRPtr->count > 0 && uRPtr < distr_end; ++uRPtr )
  {
      local_sizes.push_back( uRPtr->count );
  }
  
  pattern_t pattern( local_sizes );
  Array<Point, int, pattern_t> toSort( pattern );
  
  #ifdef DEBUG
    dash::barrier();  // only needed for better IO Output
    
    __sleep();
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": ";
    for( auto it : local_sizes ){ cout << it << ", ";}
    cout << endl;
    
    if( local_sizes.size() > myid ){
      __sleep(20);
      cout << "#" << fmt( myid, FBLUE, 2 )
           << ":"
           << " vec.size: "             << fmt( local_sizes.size(), FCYN ) // careful! may not exist -> segfault chance!
           << " veCount: "              << fmt( local_sizes[myid] , FCYN ) // careful! may not exist -> segfault chance!
           << " lsize: "                << fmt( toSort.lsize()    , FCYN )
           << " toSort lend - lbegin: " << fmt( toSort.lend( ) - toSort.lbegin( ), FCYN )
           << " Pattern size: "         << fmt( pattern.size( )   , FCYN )
           << " toSort size: "          << fmt( toSort.size  ( )  , FCYN )
           << " sizeof(Point): "        << fmt( sizeof(Point)     , FCYN )
           << " alignment: "            << fmt( (reinterpret_cast<size_t>(toSort.lbegin( )) % 64), FCYN )
           << endl;
    }
    
    #if 0
      dash::barrier();  // only needed for better IO Output
      
      if( local_sizes.size() > myid ){
        __sleep();
        uint count = 0;
        cout << "#" << myid << ": ";
        
        for( Point * i = toSort.lbegin( ); i < toSort.lend( ) ; ++i ) {
          cout << ++count << ", ";
          i->value = myid;  // test member access
          i->row   = myid;  // test member access
          i->col   = myid;  // test member access
        }
        cout << endl;
      }
    #endif
  #endif
  
  
 /* Each unit creates buckets that are send to the responsible unit later.
  * Note: every unit will get data from other units (e.g. if distribution specifies range 0-0)
  * This is a array of vectors like "vector<Point> buckets[ involvedUnits ];"
  * But on the Heap!
  */
  const size_t involvedUnits = local_sizes.size( );
  
  vector<Point> ** buckets = static_cast<vector<Point> **>(  std::malloc( involvedUnits * sizeof(vector<Point>*) )  );
  vector<Point> ** const buckets_end = buckets + involvedUnits;

  // create a vector for each bucket pointer
  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket ){ *bucket = new vector<Point>; }


  // iterate over pointsLocal and fill the buckets
  for( Point const * lclPt = pointsLocal.data(); lclPt < pointsLocal.data() + pointsLocal.size(); ++lclPt )
  {
    unitRange * uRPtr = distr;
    for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++uRPtr, ++bucket )
    {
      if( lclPt->value <= uRPtr->end )  // if value is in Range for this unit
      {
        (*bucket)->push_back(*lclPt);
        break;
      }
    }
  }
  
  #ifdef DEBUG
  dash::barrier();  // only needed for better IO Output
  
  __sleep();
  uint counter = 0;
  
  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket, ++counter )
  {
    cout << "#" << fmt( myid, FBLUE, 2 ) << ": bucket: " << fmt( counter, BBLUE ) << " :";
    
    for( auto it : **bucket ){ cout << it << ", ";}
    cout << endl;
  }
  #endif
  
 /* In comTable will be the information which unit has how many elements of type Point
  * for other units and itself.
  * The table can be thought of like "comTable[nUnits][local_sizes.size()]"
  * First every unit saves the information in row[0]
  */
  uint * comTable = static_cast<uint*>(  std::malloc( sizeof(uint) * nUnits * involvedUnits )  );
  uint * const comTable_end = comTable + ( nUnits * involvedUnits );
  
  uint * thisUnitsRow_begin     = new uint[involvedUnits];
  uint * const thisUnitsRow_end = thisUnitsRow_begin + involvedUnits;
  
  uint * poiCount = thisUnitsRow_begin;
  for( vector<Point> ** bucket = buckets; bucket < buckets_end; ++bucket, ++poiCount )
  {
    *poiCount = (*bucket)->size();
  }
  
  #ifdef DEBUG
  dash::barrier();  // only needed for better IO Output
  
  __sleep();
  
  cout << "#" << fmt( myid, FBLUE, 2 ) << ": poiCount: ";
  for( poiCount = thisUnitsRow_begin; poiCount < thisUnitsRow_end; ++poiCount )
  {
    cout << fmt( *poiCount, FRED ) << ", ";
  } cout << endl;
  #endif
  
  ret = dart_allgather(
           thisUnitsRow_begin,  // sendbuf
           comTable          ,  // recvbuf
           involvedUnits     ,  // nelem
           DART_TYPE_UINT    ,  // dtype
           team.dart_id( )   ); // team

  
  if( DART_OK != ret ) cout << "An Error while allgather has occured!" << endl; 
  
  #ifdef DEBUG
  dash::barrier();  // only needed for better IO Output
  
  __sleep();
  counter = 0;
  cout << "#" << fmt( myid, FBLUE, 2 ) << ": comTable:\n";
  for( poiCount = comTable; poiCount < comTable_end; ++poiCount, ++counter )
  {
    if( counter == involvedUnits ) counter = 0, cout << "\n";
    cout << fmt( *poiCount, FMAG ) << ", ";
  } cout << endl;
  #endif
  
  
}



int main( int argc, char* argv[] )
{  
  dash::init( &argc,&argv );
  
  myid = dash::myid( );

  ReadRowsNCols( );
  
  NArray< MATRIX_T, 2 > randMat    ( in.nrows, in.ncols );
  NArray< bool    , 2 > threshMask ( in.nrows, in.ncols );

  #ifdef DEBUG  // print error message if mask's and matrix's local size aren't identical
    if( threshMask.local_size() != randMat.local_size() )
    {
      cout << "On unit " << myid
           << " the local sizes of matrix and mask differ!\naborted on this unit\n";
      return -1;
    }
  #endif
  
  ReadMatricesAndNelts( randMat, threshMask );
  
  vector< pair<POI_T, POI_T> > result_points(nelts);
  
  winnow( in.nrows, in.ncols, randMat, threshMask, nelts, result_points );

  dash::finalize( );
}


