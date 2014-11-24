// Author: Brett Bowman

#include <map>
#include <zlib.h>
#include <stdio.h>

#include <seqan/arg_parse.h>
#include <seqan/seeds.h>
#include <seqan/index.h>

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

    addOption(parser, ArgParseOption(
                "s", "seedSize", "Size to use when searching for alignment seeds.",
                ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "seedSize", "12");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Extract the arguments from the Parser
    std::string query;
    std::string reference;
    int seedSize;

    getArgumentValue(query, parser, 0);
    getArgumentValue(reference, parser, 1);
    getOptionValue(seedSize, parser, "seedSize");

    // Use the options to set the configs and scoring schemes
    Score<long, Simple> scoringScheme(4, -13, -7);

    typedef FindSeedsConfig<12> TConfig;

    // Read the reference sequences into memory
    ReferenceSet refSet = ReferenceSet(reference);
    auto refSetIndex = refSet.GetIndex<TConfig>();
    indexRequire(refSetIndex, QGramSADir());  // On-demand index creation.
   
    //typedef StringSet<Dna5String> TStringSet;
    //typedef Index<TStringSet, TConfig::IndexType> TIndex;
    //typedef Iterator<TIndex, TopDown<ParentLinks<>>> TIter;

    //TIter::Type it(refSetIndex);
    //do {
    //    // Print theletters from the root to the current node
    //    std::cout << representative(it) << std::endl;

    //    if (!goDown(it) && !goRight(it))
    //        while (goUp(it) && !goRight(it));
    //} while (!isRoot(it));

    //TShape & shape = indexShape(index);
    //hashInit(shape, begin(query, Standard()));
    //for (TIterator it = begin(query, Standard()); it != end(query, Standard()) - 12; ++it)
    //{
    //    std::cout << "Occ at: ";
    //    hashNext(shape, it);
    //    for (unsigned i = 0; i < length(getOccurrences(index, shape)); ++i)
    //        std::cout << getOccurrences(index, shape)[i] << " ";
    //    std::cout << std::endl;
    //}

    // Create an iterator for the query sequences and a pair for it to return to
    SequenceReader seqReader = SequenceReader(query);
    std::pair<size_t, SequenceRecord> idxAndRecord;

    // Define the variable where the initial hits for the query will be stored
    //map<size_t, SeedSet<Simple>> querySeedHits;
    SeedSet<Simple> querySeedHits;

    for ( ; seqReader.GetNext(idxAndRecord) ; )
    {
        std::cout << "Query #" << idxAndRecord.first+1 << " - "<< idxAndRecord.second.Id << std::endl;

        // For each query sequence, compare it to the ReferenceSet
        //for (unsigned i = 0; i < refSet.Length(); ++i)
        //    std::cout << "Ref " << refSet.Ids()[i] << '\t' << refSet.Sequences()[i] << std::endl;

        // Find the Kmer matches for the current query sequence
        FindSeeds3(querySeedHits, refSetIndex, idxAndRecord.second.Seq);
        std::cout << "Finished finding seeds" << std::endl;

        // Chain the initial Kmer hits into an alignment
        auto align = SeedsToAlignment(idxAndRecord.second.Seq, 
                                      refSet.Sequences()[0],
                                      querySeedHits,
                                      scoringScheme);

        //querySeedHits.clear();
    }
}
