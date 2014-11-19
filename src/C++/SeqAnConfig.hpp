
#pragma once

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

// Config for finding alignment seeds with suffix arrays / QGrams 
//
// TODO (lhepler) : investigate default values other than 10
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
typedef AlignConfig<false, false, false, false> GlobalAlignConfig();
