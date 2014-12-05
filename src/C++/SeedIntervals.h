
#pragma once

#include <vector>

#include <seqan/seeds.h>

#include "SeqAnConfig.hpp"

using namespace std;
using namespace seqan;

// These tuples represent a SeedSet and a Start and End index that
//   identify the seeds representing the 5' and 3' most seed anchors
//   for a hit
typedef std::tuple<size_t, size_t, size_t> SeedInterval;

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
