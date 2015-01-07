#include <vector>

#include <seqan/seeds.h>

#include "utils/Seed.cpp"
#include "utils/SeedChain.cpp"
#include "utils/ReferencedSeedChain.cpp"
#include "SeedIntervals.hpp"
#include "ReferenceSet.hpp"
#include "config/SeqAnConfig.hpp"

using namespace seqan;
using namespace srsli;

template<typename TSeed>
void SeedSetsToSeedVectors(std::vector<TSeedSet>& sets,
                           std::vector<std::vector<TSeed>>& vectors)
{
    typedef seqan::Iterator<TSeedSet>::Type TIter;

    // Iterate over each set, each containing the hits to one reference
    for (size_t i = 0; i < sets.size(); ++i)
    {
        // Initialize destination Seed vector to the appropriate size
        vectors[i] = std::vector<TSeed>( length(sets[i]) );

        // Iterate over each seed, adding it to the vector
        size_t j = 0;
        for (TIter it = begin(sets[i], seqan::Standard()); it != end(sets[i], seqan::Standard()); ++it)
        {
            vectors[i][j] = *it;
            ++j;
        }

        // Sort each vector after it's been filled by reference position
        std::sort(vectors[i].begin(), vectors[i].end(), seedReferencePositionComparer);
    }
}


template<typename TSeed>
int AdvanceIndexToIntervalEnd(const std::vector<TSeed>& seeds,
                              const size_t& nSeeds,
                              const size_t& maxIntervalSize,
                              const size_t& start,
                              size_t& end)
{
    // Calculate the maximal interval end position once
    size_t maxEndPos = beginPositionV(seeds[start]) + maxIntervalSize;

    while (// While we are not at the end of the vector
           end+1 < nSeeds and
           // ... and the current seed ends before the maximal end position
           maxEndPos > endPositionV(seeds[end+1])) {

        // increment the index of the end seed by one
        end++;
    }

    // If we made it this far, return 0 for successful completion
    return 0;
}


template<typename TSeed>
int GetSeedIntervals(std::vector<SeedInterval>& intervals,
                     const std::vector<std::vector<TSeed>>& seedHits,
                     const size_t& maxIntervalSize)
{
    // Iterate over each SeedSet, looking for intervals in each
    for (size_t refIdx = 0; refIdx < seedHits.size(); ++refIdx)
    {
        size_t nSeeds = length(seedHits[refIdx]);
        size_t currEndIdx = 1, prevEndIdx = 0;

        // Iterate over each seed, treating it as a possible interval start index
        for (size_t startIdx = 0; startIdx < nSeeds; ++startIdx)
        {
            // Find the last possible Seed for an interval starting at startIdx
            AdvanceIndexToIntervalEnd(seedHits[refIdx], nSeeds, maxIntervalSize, startIdx, currEndIdx);

            // If the end index hasn't move skip to the next iteration
            if (currEndIdx == prevEndIdx)
                continue;

            // Otherwise save the current interval and terminal index
            SeedInterval currInterval(refIdx, startIdx, currEndIdx);
            intervals.push_back(currInterval);
            prevEndIdx = currEndIdx;

            // If the current interval ends with the last seed, there
            //    are no more possible new intervals to find, so 
            //    we quit
            if (currEndIdx == nSeeds-1)
                break;
        }
    }

    // If we made it this far, return 0 for successful completion
    return 0;
}


template<typename TSeed>
int SeedSetFromSeedInterval(TSeedSet& seedSet,
                            const SeedInterval& interval,
                            const std::vector<TSeed>& seeds)
{
    size_t start = std::get<1>(interval);
    size_t end = std::get<2>(interval);
    for (size_t i = start; i < end; ++i)
    {
        // Note that since these seeds were merged when they were
        //    originally found, we can add them simply this time
        addSeed(seedSet, seeds[i], Single());
    }

    // If we made it this far, return 0 for successful completion
    return 0;
}


template<typename TSeed>
int SeedIntervalsToSeedChains(std::vector<ReferencedSeedChain>& chains,
                              const std::vector<std::vector<TSeed>>& seedVecs,
                              const std::vector<SeedInterval>& intervals)
{
    size_t minSeedChainBases = 30;

    // Allocate a seedSet and chain for intermediate use
    TSeedSet seedSet;
    TSeedChain chain; 
    ReferencedSeedChain refChain;
    size_t prevIdx = 0;
    int endPos, prevEndPos;
    int startPos, prevStartPos;
    size_t numChainBases;

    for (size_t i = 0; i < intervals.size(); ++i)
    {
        // First identify the Reference / Anchor list we will be using
        ReferencedSeedChain refChain;
        refChain.referenceIndex = std::get<0>(intervals[i]);

        // Fill the seed set with the appropriate seeds and Chain
        clear(seedSet);
        SeedSetFromSeedInterval(seedSet, 
                                intervals[i], 
                                seedVecs[refChain.referenceIndex]);

        // Chain the seeds together and find the chain's start and end points
        chainSeedsGlobally(refChain.chain, seedSet, SparseChaining());
     
        // Skip seed chains with very little supporting evidence
        if (minSeedChainBases > SumSeedChainBases(refChain.chain))
            continue;

        // If we pass that first filter, we calculate Start and End positions
        startPos = GetSeedChainStartPos(chain);
        endPos = GetSeedChainEndPos(chain);

        // If this is the first chain
        if (chains.size() == 0)
        {
            // .. we automatically keep the chain and its Start and End positions
            chains.push_back(refChain);
            prevStartPos = startPos;
            prevEndPos = endPos;

        // If we have the same start but a longer end, replace the previous chain
        } else if (startPos == prevStartPos && endPos > prevEndPos ) {
            chains[prevIdx] = refChain;

        // On the other hand, if we have the same end and a lesser-or-matching start, skip it
        } else if (startPos >= prevStartPos && endPos == prevEndPos) {
            continue;

        // Finally, if neither of these is true and we have a novel chain, add it
        } else {
            chains.push_back(refChain);
            prevStartPos = startPos;
            prevEndPos = endPos;
            prevIdx++;
        }
    }

    // Finally, we sort the SeedChains we found by the number of bp matches they represent
    std::sort(chains.begin(), chains.end(), refSeedChainNumBasesComparer);

    // If we made it this far, return 0 for successful completion
    return 0;
}
