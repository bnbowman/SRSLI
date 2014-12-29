#include <vector>

#include <seqan/seeds.h>

#include "utils/Seed.cpp"
#include "utils/SeedChain.cpp"
#include "SeedIntervals.hpp"
#include "ReferenceSet.hpp"
#include "config/SeqAnConfig.hpp"

using namespace seqan;
using namespace srsli;

template<typename TSeed>
void SeedSetsToSeedVectors(vector<TSeedSet>& sets,
                           vector<vector<TSeed>>& vectors)
{
    typedef seqan::Iterator<TSeedSet>::Type TIter;

    // Iterate over each set, each containing the hits to one reference
    for (size_t i = 0; i < sets.size(); ++i)
    {
        // Initialize destination Seed vector to the appropriate size
        vectors[i] = vector<TSeed>( length(sets[i]) );

        // Iterate over each seed, adding it to the vector
        size_t j = 0;
        for (TIter it = begin(sets[i], seqan::Standard()); it != end(sets[i], seqan::Standard()); ++it)
        {
            vectors[i][j] = *it;
            ++j;
        }

        // Sort each vector after it's been filled by reference position
        std::sort(vectors[i].begin(), vectors[i].end(), seedReferencePositionComparer);

        //for (int j = 0; j < vectors[i].size(); ++j)
        //{
        //    std::cout << vectors[i][j] << std::endl;
        //}
    }
}


template<typename TSeed>
int AdvanceIndexToIntervalEnd(const vector<TSeed>& seeds,
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
int GetSeedIntervals(vector<SeedInterval>& intervals,
                     const vector<vector<TSeed>>& seedSets,
                     const size_t& maxIntervalSize)
{
    // Iterate over each SeedSet, looking for intervals in each
    for (size_t i = 0; i < seedSets.size(); ++i)
    {
        size_t nSeeds = length(seedSets[i]);
        size_t currEndIdx = 1, prevEndIdx = 0;

        // Iterate over each seed, treating it as a possible interval start index
        for (size_t startIdx = 0; startIdx < nSeeds; ++startIdx)
        {
            // Find the last possible Seed for an interval starting at startIdx
            AdvanceIndexToIntervalEnd(seedSets[i], nSeeds, maxIntervalSize, startIdx, currEndIdx);

            // If the end index hasn't move skip to the next iteration
            if (currEndIdx == prevEndIdx)
                continue;

            // Otherwise save the current interval and terminal index
            SeedInterval currInterval(i, startIdx, currEndIdx);
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
                            const vector<TSeed>& seeds)
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
int SeedIntervalsToSeedChains(vector<vector<TSeed>>& chains,
                              const vector<vector<TSeed>>& seedVecs,
                              const vector<SeedInterval>& intervals)
{
    size_t minSeedChainBases = 30;

    // Allocate a seedSet and chain for intermediate use
    TSeedSet seedSet;
    vector<TSeed> chain;
    size_t prevIdx;
    int endPos, prevEndPos;
    int startPos, prevStartPos;
    size_t numChainBases;

    for (size_t i = 0; i < intervals.size(); ++i)
    {
        // First identify the Reference / Anchor list we will be using
        size_t seedVecIdx = std::get<0>(intervals[i]);

        // Fill the seed set with the appropriate seeds and Chain
        clear(seedSet);
        SeedSetFromSeedInterval(seedSet, intervals[i], seedVecs[seedVecIdx]);

        // Chain the seeds together and find the chain's start and end points
        chain.clear();
        chainSeedsGlobally(chain, seedSet, SparseChaining());
        startPos = GetSeedChainStartPos(chain);
        endPos = GetSeedChainEndPos(chain);
     
        // Skip seed chains with very little supporting evidence
        if (minSeedChainBases > SumSeedChainBases(chain))
            continue;

        // Add the chain to our current vector of possible chains
        //    which is automatic if there are no previous chains
        if (chains.size() == 0)
        {
            chains.push_back(chain);
            prevStartPos = startPos;
            prevEndPos = endPos;
            prevIdx = 0;

        // If we have the same start but a longer end, replace the previous chain
        } else if (startPos == prevStartPos && endPos > prevEndPos ) {
            chains[prevIdx] = chain;

        // On the other hand, if we have the same end and a lesser-or-matching start, skip it
        } else if (startPos >= prevStartPos && endPos == prevEndPos) {
            continue;

        // Finally, if neither of these is true and we have a novel chain, add it
        } else {
            chains.push_back(chain);
            prevStartPos = startPos;
            prevEndPos = endPos;
            prevIdx++;
        }
    }

    // Finally, we sort the SeedChains we found by the number of bp matches they represent
    std::sort(chains.begin(), chains.end(), seedChainNumBasesComparer);

    // If we made it this far, return 0 for successful completion
    return 0;
}
