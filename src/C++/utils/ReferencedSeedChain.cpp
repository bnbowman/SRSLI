#pragma once
#include <vector>
#include <seqan/seeds.h>

#include "../config/SeqAnConfig.hpp"
#include "../config/Types.hpp"
#include "Seed.cpp"

// Utility Functions
size_t SumReferencedSeedChainBases(const ReferencedSeedChain& refChain)
{
    const TSeedChain* chain = &refChain.chain;
    size_t sum = 0;
    for (size_t i = 0; i < length(*chain); ++i)
    {
        sum += GetSeedNumBases(chain->at(i));
    }
    return sum;
}

// Functors
struct ReferencedSeedChainNumBasesFunctor {
    bool operator()(const ReferencedSeedChain& chain1, const ReferencedSeedChain& chain2)
    {
        return SumReferencedSeedChainBases(chain1) > SumReferencedSeedChainBases(chain2);
    }
} refSeedChainNumBasesComparer;
// End Functors
