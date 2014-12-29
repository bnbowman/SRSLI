#pragma once

#include <vector>

#include "SeqAnConfig.hpp"

using namespace seqan;

// These tuples represent a SeedSet and a Start and End index that
//   identify the seeds representing the 5' and 3' most seed anchors
//   for a hit
typedef std::tuple<size_t, size_t, size_t> SeedInterval;

// This pair represents an ordered vector of seeds and the index
//   of the reference sequence from which it came
typedef std::pair<size_t, TSeedChain> ReferencedSeedChain;
