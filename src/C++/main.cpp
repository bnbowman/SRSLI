// Author: Brett Bowman

#include <math.h>
#include <zlib.h>
#include <stdio.h>

#include <seqan/arg_parse.h>
#include <seqan/seeds.h>
#include <seqan/index.h>

#include "Utils.hpp"
#include "Version.hpp"
#include "SeqAnConfig.hpp"
#include "SrsliParameters.h"
#include "ReferenceSet.hpp"
#include "SequenceReader.hpp"
#include "SparseAlignment.hpp"
#include "SeedIntervals.cpp"

using namespace seqan;
using namespace srsli;

// Entry point
int main(int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("srsli");
    setDate(parser, Version::Date());
    setVersion(parser, Version::VersionString());
    setShortDescription(parser, "Successively Refined alignment after Sparse Local Indexing");
    addDescription(parser, "Quickly and efficiently map single-molecule sequencing reads to a "
                           "set of reference sequences via sparse local indexing, then refine the "
                           "results with standard dynamic programming methods");
    addUsageLine(parser,  "\\fIQUERY\\fP \\fIREFERENCE\\fP [\\fIOPTIONS\\fP]");

    // Define Required (Positional) Arguments
    addArgument(parser, ArgParseArgument(
                ArgParseArgument::STRING, "QUERY"));
    addArgument(parser, ArgParseArgument(
                ArgParseArgument::STRING, "REFERENCE"));

    // Define Optional arguments
    addOption(parser, ArgParseOption(
                "m", "minScore", "Minimum score to require from any alignment",
                ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption(
                "n", "nCandidates", "Number of candidate alignments to score",
                ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption(
                "s", "seedSize", "Size to use when searching for alignment seeds.",
                ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption(
                "v", "verbosity", "Verbosity of process information to report [0..3].",
                ArgParseArgument::INTEGER, "INT"));

    // Set default values
    setDefaultValue(parser, "minScore",    "1000");
    setDefaultValue(parser, "nCandidates", "1"   );
    setDefaultValue(parser, "seedSize",    "12"  );
    setDefaultValue(parser, "verbosity",   "1"   );

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // If parsing was successful, initialize the parameters
    SrsliParameters params;

    // Extract the arguments from the Parser
    std::string query;
    std::string reference;
    int minScore;
    int nCandidates;
    int seedSize;
    int verbosity;

    getArgumentValue(query,     parser, 0);
    getArgumentValue(reference, parser, 1);
    getOptionValue(minScore,    parser, "minScore");
    getOptionValue(nCandidates, parser, "nCandidates");
    getOptionValue(seedSize,    parser, "seedSize");
    getOptionValue(verbosity,   parser, "verbosity");

    // Use the options to set the configs and scoring schemes
    Score<long, Simple> scoringScheme(4, -13, -7);

    typedef FindSeedsConfig<12> TConfig;

    // Read the reference sequences into memory
    ReferenceSet refSet = ReferenceSet(reference);
    auto refSetIndex = refSet.GetIndex<TConfig>();
    indexRequire(refSetIndex, QGramSADir());  // On-demand index creation.
   
    // Create an iterator for the query sequences and a pair for it to return to
    SequenceReader seqReader = SequenceReader(query);
    std::pair<size_t, SequenceRecord> idxAndRecord;

    // Define the variable where the initial hits for the query will be stored
    vector<TSeedSet> querySeedSets(2, TSeedSet());
    vector<vector<TSeed>> querySeedHits(2);
    vector<vector<TSeed>> querySeedChains;
    //TSeedSet querySeedHits;

    std::cout << refSet.Size() << std::endl;
    std::cout << refSet.Length() << std::endl;

    for ( ; seqReader.GetNext(idxAndRecord) ; )
    {
        std::cout << "Query #" << idxAndRecord.first+1 << " - "<< idxAndRecord.second.Id << std::endl;

        // Calculate the maximum expected interval size, given the query length
        size_t maxIntervalLength = length(idxAndRecord.second.Seq) * params.maxNetIndelRate;

        // Find the Kmer matches for the current query sequence
        std::cout << "Looking for seeds" << std::endl;
        FindSeeds3(querySeedSets, refSetIndex, refSet.Size(), idxAndRecord.second.Seq);
        std::cout << "Finished finding seeds" << std::endl;

        // Sort the Kmer matches by reference position, which involves converting them to vectors
        std::cout << "Sorting seeds" << std::endl;
        SeedSetsToSeedVectors(querySeedSets, querySeedHits);
        std::cout << "Finished sorting seeds" << std::endl;

        std::cout << "Selecting seed intervals" << std::endl;
        vector<SeedInterval> seedIntervals;
        GetSeedIntervals(seedIntervals, querySeedHits, maxIntervalLength);
        std::cout << "Finished selecting seed intervals" << std::endl;

        std::cout << "Converting seed intervals to seed chains" << std::endl;
        SeedIntervalsToSeedChains(querySeedChains, querySeedHits, seedIntervals);
        std::cout << "Finished converting seed intervals to chains" << std::endl;

        for (size_t i = 0; i < querySeedChains.size(); ++i)
        {
            auto temp = querySeedChains[i];
            std::cout << "Idx #" << i << " L " << temp.size() << " Pos "<< GetSeedChainStartPos(temp) << " -- " << GetSeedChainEndPos(temp) << std::endl;
        }

        //vector<int> seedSetOrder = GetSeedSetOrder(querySeedHits);

        //int maxAligns = std::min((int)seedSetOrder.size(), nCandidates);
        //for (size_t i = 0; i < maxAligns; ++i)
        //{
            // Chain the initial Kmer hits into an alignment
        //    int seedSetIdx = seedSetOrder[i];
        //    auto align = SeedsToAlignment(idxAndRecord.second.Seq, 
        //                                  refSet.Sequences()[seedSetIdx],
        //                                  querySeedHits[seedSetIdx],
        //                                  scoringScheme);
        
        //    std::cout << align << std::endl;
        //}

        // Empty the seedHits variable before the next iteration
        querySeedHits.clear();
    }
}
