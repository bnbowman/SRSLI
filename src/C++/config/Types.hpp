#pragma once

#include <vector>

#include "SeqAnConfig.hpp"

using namespace seqan;

// These tuples represent a SeedSet and a Start and End index that
//   identify the seeds representing the 5' and 3' most seed anchors
//   for a hit
typedef std::tuple<size_t, size_t, size_t> SeedInterval;

// This struct represents an ordered vector of seeds and the index
//   of the reference sequence from which it came
struct ReferencedSeedChain {
    size_t referenceIndex;
    TSeedChain chain;

    ReferencedSeedChain() {}

    ReferencedSeedChain(size_t i, TSeedChain c)
        : referenceIndex( i )
        , chain( c )
    {}
};

// This struct represents the Start and End positions within a
//   pair of sequences to be aligned
struct region_t {
    int queryStart;
    int queryEnd;
    int refStart;
    int refEnd;
};

struct ReferenceRecord {
    CharString id;
    TDna* seq;
    int orientation;

    ReferenceRecord() {}

    ReferenceRecord(CharString i, TDna* s, int o)
        : id( i )
        , seq( s )
        , orientation( o )
    {}
};
