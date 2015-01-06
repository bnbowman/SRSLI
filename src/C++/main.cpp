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

#include "config/SeqAnConfig.hpp"
#include "parameters/SrsliParameters.hpp"
#include "Headers.cpp"
#include "ReferenceSet.hpp"
#include "AlignmentRecord.hpp"
#include "SequenceReader.hpp"
#include "FindSeeds.hpp"
#include "SeedIntervals.cpp"
#include "SparseAlignment.cpp"

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
    std::vector<TSeedSet> querySeedSets(2, TSeedSet());
    std::vector<std::vector<TSeed>> querySeedHits(2);
    std::vector<ReferencedSeedChain> querySeedChains;

    // Initialize the vector that will get the alignment results
    std::vector<AlignmentRecord> results;

    for ( ; seqReader.GetNext(idxAndRecord) ; ) {
        // Display current query
        SequenceRecord* record = &idxAndRecord.second;
        std::cout << "Query #" << idxAndRecord.first+1
                  << " - " << record->Id << std::endl;

        // Calculate the maximum expected interval size, given the query length
        size_t maxIntervalLength = length(record->Seq) * params.maxNetIndelRate;

        // Find the Kmer matches for the current query sequence
        std::cout << "Looking for seeds" << std::endl;
        FindSeeds(querySeedSets, refSetIndex, refSet.Size(), record->Seq);
        std::cout << "Finished finding seeds" << std::endl;

        // Sort the Kmer matches by reference position, which involves
        //    converting them to vectors
        std::cout << "Sorting seeds" << std::endl;
        SeedSetsToSeedVectors(querySeedSets, querySeedHits);
        std::cout << "Finished sorting seeds" << std::endl;

        std::cout << "Selecting seed intervals" << std::endl;
        std::vector<SeedInterval> seedIntervals;
        GetSeedIntervals(seedIntervals, querySeedHits, maxIntervalLength);
        std::cout << "Finished selecting seed intervals" << std::endl;

        std::cout << "Converting seed intervals to seed chains" << std::endl;
        SeedIntervalsToSeedChains(querySeedChains,
                                  querySeedHits,
                                  seedIntervals);
        std::cout << "Finished conversion to seed chains" << std::endl;

        int maxAligns = std::min((int)querySeedChains.size(), params.nCandidates);

        // Chain the initial Kmer hits into an alignment
        RefChainsToAlignments(results,
                              record->Seq,
                              refSet,
                              querySeedChains,
                              scoringScheme,
                              maxAligns,
                              params.minAccuracy,
                              params.maxChainBuffer,
                              params.alignmentAnchor);

        // Empty the seedHits variable before the next iteration
        querySeedHits.clear();
    }

    std::cout << results.size() << std::endl;
    std::cout << M1Header << std::endl;
    for (size_t i = 0; i < results.size(); ++i)
    {
        std::cout << results[i].toM1Record() << std::endl;
    }
}
