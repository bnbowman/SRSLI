// Copyright 2014 Brett Bowman

#include <math.h>
#include <zlib.h>
#include <stdio.h>
#include <utility>
#include <string>
#include <vector>

#include <seqan/arg_parse.h>
#include <seqan/seeds.h>
#include <seqan/index.h>

#include "Utils.hpp"
#include "Version.hpp"
#include "SeqAnConfig.hpp"
#include "SrsliParameters.hpp"
#include "ReferenceSet.hpp"
#include "SequenceReader.hpp"
#include "SparseAlignment.hpp"
#include "SeedIntervals.cpp"

using namespace seqan;
using namespace srsli;

// Entry point
int main(int argc, char const ** argv) {

    // Parse the command-line arguments into usable parameters
    SrsliParameters params(argc, argv);

    // Abort if we could not parse the supplied arguments
    if (params.parseOk == 0)
        return 1;

    // Use the options to set the configs and scoring schemes
    Score<int64_t, Simple> scoringScheme(4, -13, -7);

    // Read the reference sequences into memory
    typedef FindSeedsConfig<12> TConfig;
    ReferenceSet refSet = ReferenceSet(params.reference);
    auto refSetIndex = refSet.GetIndex<TConfig>();
    indexRequire(refSetIndex, QGramSADir());  // On-demand index creation.

    // Create an iterator for the query sequences and a pair for it to return to
    SequenceReader seqReader = SequenceReader(params.query);
    std::pair<size_t, SequenceRecord> idxAndRecord;
    // Define the variable where the initial hits for the query will be stored
    vector<TSeedSet> querySeedSets(2, TSeedSet());
    vector<vector<TSeed>> querySeedHits(2);
    vector<vector<TSeed>> querySeedChains;

    for ( ; seqReader.GetNext(idxAndRecord) ; ) {
        // Display current query
        SequenceRecord* record = &idxAndRecord.second;
        std::cout << "Query #" << idxAndRecord.first+1
                  << " - " << record->Id << std::endl;

        // Calculate the maximum expected interval size, given the query length
        size_t maxIntervalLength = length(record->Seq) * params.maxNetIndelRate;

        // Find the Kmer matches for the current query sequence
        std::cout << "Looking for seeds" << std::endl;
        FindSeeds3(querySeedSets, refSetIndex, refSet.Size(), record->Seq);
        std::cout << "Finished finding seeds" << std::endl;

        // Sort the Kmer matches by reference position, which involves
        //    converting them to vectors
        std::cout << "Sorting seeds" << std::endl;
        SeedSetsToSeedVectors(querySeedSets, querySeedHits);
        std::cout << "Finished sorting seeds" << std::endl;

        std::cout << "Selecting seed intervals" << std::endl;
        vector<SeedInterval> seedIntervals;
        GetSeedIntervals(seedIntervals, querySeedHits, maxIntervalLength);
        std::cout << "Finished selecting seed intervals" << std::endl;

        std::cout << "Converting seed intervals to seed chains" << std::endl;
        SeedIntervalsToSeedChains(querySeedChains,
                                  querySeedHits,
                                  seedIntervals);
        std::cout << "Finished conversion to seed chains" << std::endl;

        for (size_t i = 0; i < querySeedChains.size(); ++i) {
            auto temp = querySeedChains[i];
            std::cout << "Idx #" << i
                      << " L " << temp.size()
                      << " Pos "<< GetSeedChainStartPos(temp)
                      << " -- " << GetSeedChainEndPos(temp) << std::endl;
        }

        // vector<int> seedSetOrder = GetSeedSetOrder(querySeedHits);
        // int maxAligns = std::min((int)seedSetOrder.size(), nCandidates);
        // for (size_t i = 0; i < maxAligns; ++i)
        // {
        // Chain the initial Kmer hits into an alignment
        //    int seedSetIdx = seedSetOrder[i];
        //    auto align = SeedsToAlignment(idxAndRecord.second.Seq,
        //                                  refSet.Sequences()[seedSetIdx],
        //                                  querySeedHits[seedSetIdx],
        //                                  scoringScheme);
        //    std::cout << align << std::endl;
        // }
        // Empty the seedHits variable before the next iteration
        querySeedHits.clear();
    }
}
