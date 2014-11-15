// Author: Brett Bowman

#include <map>
#include <zlib.h>
#include <stdio.h>

#include <seqan/arg_parse.h>
#include <seqan/seeds.h>

#include "Utils.hpp"
#include "Version.hpp"
#include "SeqAnConfig.hpp"
#include "ReferenceSet.hpp"
#include "SequenceReader.hpp"
#include "SparseAlignment.hpp"

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

    addArgument(parser, ArgParseArgument(
                ArgParseArgument::STRING, "QUERY"));
    addArgument(parser, ArgParseArgument(
                ArgParseArgument::STRING, "REFERENCE"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract the arguments from the Parser
    std::string query;
    std::string reference;

    getArgumentValue(query, parser, 0);
    getArgumentValue(reference, parser, 1);

    // Read the reference sequences into memory
    ReferenceSet refSet = ReferenceSet(reference);
    auto refSetIndex = refSet.GetIndex<FindSeedsConfig<12>>();

    // Create an iterator for the query sequences and a pair for it to return to
    SequenceReader seqReader = SequenceReader(query);
        std::pair<size_t, SequenceRecord> idxAndRecord;

    // Define the variable where the initial hits for the query will be stored
    map<size_t, SeedSet<Simple>> querySeedHits;

    for ( ; seqReader.GetNext(idxAndRecord) ; )
    {
        std::cout << "Query " << idxAndRecord.first << std::endl;

        // For each query sequence, compare it to the ReferenceSet
        //for (unsigned i = 0; i < refSet.Length(); ++i)
        //    std::cout << "Ref " << refSet.Ids()[i] << '\t' << refSet.Sequences()[i] << std::endl;

        FindSeeds(querySeedHits, refSetIndex, idxAndRecord.second.Seq);

        querySeedHits.clear();
    }
}
