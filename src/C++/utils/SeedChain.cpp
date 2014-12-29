#pragma once
#include <vector>
#include <seqan/seeds.h>

#include "../SeqAnConfig.hpp"
#include "Seed.cpp"


// Utility Functions
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

template<typename TSeedChain>
int GetSeedChainStartPos(const TSeedChain& chain)
{
    return beginPositionV( chain[0] );
}

template<typename TSeedChain>
int GetSeedChainEndPos(const TSeedChain& chain)
{
    return endPositionV( chain[chain.size()-1] );
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
