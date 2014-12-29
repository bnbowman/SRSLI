
#pragma once

#include <vector>

#include <seqan/seeds.h>

#include "config/SeqAnConfig.hpp"
#include "config/Types.hpp"

using namespace std;
using namespace seqan;

template<typename TSeed>
void SeedSetsToSeedVectors(vector<TSeedSet>& sets,
        vector<vector<TSeed>>& strings);

template<typename TSeed>
int AdvanceIndexToIntervalEnd(const TSeedSet& seedSet,
        const size_t& nSeeds,
        const size_t& maxIntervalSize,
        const size_t& start,
        size_t& end);

template<typename TSeed>
int  GetSeedIntervals(vector<SeedInterval>& intervals,
        vector<TSeedSet>& seedSets);
