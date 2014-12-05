
#pragma once

#include <vector>

#include <seqan/seeds.h>

#include "SeqAnConfig.hpp"

using namespace std;
using namespace seqan;

typedef std::tuple<size_t, size_t, size_t> SeedInterval;

template<typename T>
void SeedSetsToSeedVectors(const vector<TSeedSet>& sets,
                           vector<vector<TSeed>>& strings)

template<typename T>
int AdvanceIndexToIntervalEnd(const TSeedSet& seedSet,
                              const size_t& maxIntervalSize,
                              const size_t& start,
                              size_t& end)

template<typename T>
int  GetSeedIntervals(vector<SeedInterval>& intervals,
                      vector<TSeedSet>& seedSets)
