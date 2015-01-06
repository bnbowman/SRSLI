
#pragma once

#include <vector>

#include <seqan/seeds.h>

#include "config/SeqAnConfig.hpp"
#include "config/Types.hpp"

using namespace seqan;

template<typename TSeed>
void SeedSetsToSeedVectors(std::vector<TSeedSet>& sets,
                           std::vector<std::vector<TSeed>>& strings);

template<typename TSeed>
int AdvanceIndexToIntervalEnd(const TSeedSet& seedSet,
                              const size_t& nSeeds,
                              const size_t& maxIntervalSize,
                              const size_t& start,
                              size_t& end);

template<typename TSeed>
int  GetSeedIntervals(std::vector<SeedInterval>& intervals,
                      std::vector<TSeedSet>& seedSets);
