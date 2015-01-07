#pragma once
#include <vector>
#include <seqan/seeds.h>

#include "../config/SeqAnConfig.hpp"
#include "Seed.cpp"


// Utility Functions
template<typename TSeedChain>
int beginPositionV(const TSeedChain& chain)
{
    return beginPositionV( front(chain) );
}

template<typename TSeedChain>
int endPositionV(const TSeedChain& chain)
{
    return endPositionV( back(chain) );
}

template<typename TSeedChain>
int beginPositionH(const TSeedChain& chain)
{
    return beginPositionH( front(chain) );
}

template<typename TSeedChain>
int endPositionH(const TSeedChain& chain)
{
    return endPositionH( back(chain) );
}

template<typename TSeedChain>
size_t SumSeedChainBases(const TSeedChain& chain)
{
    size_t sum = 0;
    for (size_t i = 0; i < length(chain); ++i)
    {
        sum += GetSeedNumBases(chain[i]);
    }
    return sum;
}
// End Utility Functions

// Functors
struct SeedChainNumBasesFunctor {
    bool operator()(const TSeedChain& chain1, const TSeedChain& chain2)
    {
        return SumSeedChainBases(chain1) > SumSeedChainBases(chain2);
    }
} seedChainNumBasesComparer;
// End Functors
