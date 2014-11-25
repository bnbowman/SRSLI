
#pragma once

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

// Declare our Seed Config and associated Seed type
struct SeedConfig
{
    typedef size_t TPosition;
    //typedef size_t TScore;
    typedef MakeSigned_<size_t>::Type TDiagonal;
    typedef size_t TSize;
    typedef float TScoreValue;   // Changed to Float to be compatible with BLASR seed scoring
};

//typedef Seed<SeedConfig> TSeed;
//typedef SeedSet<SeedConfig> TSeedSet;
//typedef String<SeedConfig> TSeedString;

typedef Seed<Simple> TSeed;
typedef SeedSet<Simple> TSeedSet;
typedef String<TSeed> TSeedString;

// Config for finding alignment seeds with suffix arrays / QGrams 
//
// TODO (lhepler) : investigate default values other than 10
//template <size_t TSize = 12, typename TShape = UngappedShape<TSize>, typename TIndex = FMIndex<>>
template <size_t TSize = 12, typename TShape = UngappedShape<TSize>, typename TIndex = IndexQGram<TShape>>
class FindSeedsConfig
{
    public:
        typedef TShape ShapeType;
        typedef TIndex IndexType;
        static const int Size = TSize;
};

const FindSeedsConfig<> DefaultFindSeedsConfig();

// Alignment configurations
typedef AlignConfig<false, true, true, false> GlobalAlignConfig();
