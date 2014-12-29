#pragma once
#include <seqan/seeds.h>
#include "../SeqAnConfig.hpp"


// Begin Utility Functions
template<typename TSeed>
size_t GetSeedNumBases(const TSeed seed)
{
    return endPositionV(seed) - beginPositionV(seed);
}
// End Utility Functions


// Begin Functors
struct SeedReferencePositionFunctor {
    bool operator()(const TSeed& seed1, const TSeed& seed2)
    {
        return beginPositionV(seed1) < beginPositionV(seed2);
    }
} seedReferencePositionComparer;
// End Functors
