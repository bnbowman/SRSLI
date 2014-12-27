// Author: Brett Bowman

#pragma once

using namespace srsli;

class SrsliParameters {
    public:
        // Constructor 
        SrsliParameters(int, char const **);

        // Positional (Required) Arguments
        std::string query;
        std::string reference;

        // Optional Arguments
        int minScore;
        int nCandidates;
        int seedSize;
        int verbosity;

        // Hidden and fixed parameters
        float maxNetIndelRate;
        int parseOk;

    private:
        seqan::ArgumentParser parser;
        seqan::ArgumentParser::ParseResult result;

    void Initialize() 
    {
        if (parseOk == 0)
            return;

        // Extract Positional Arguments into variables
        getArgumentValue(query,     parser, 0);
        getArgumentValue(reference, parser, 1);

        // Extract Optional Arguments into variables
        getOptionValue(minScore,    parser, "minScore");
        getOptionValue(nCandidates, parser, "nCandidates");
        getOptionValue(seedSize,    parser, "seedSize");
        getOptionValue(verbosity,   parser, "verbosity");

        // Set hiddeen parameters
        maxNetIndelRate = 1.30;
    }

    seqan::ArgumentParser SetupParser()
    {
        // Setup ArgumentParser.
        seqan::ArgumentParser parser("srsli");
        setDate(parser, Version::Date());
        setVersion(parser, Version::VersionString());
        setShortDescription(parser,
                "Successively Refined alignment after Sparse Local Indexing");
        addDescription(parser, "Quickly and efficiently map single-molecule"
                " sequencing reads to a set of reference sequences via sparse"
                " local indexing, then refine the results with standard dynamic"
                " programming methods");
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
                "v", "verbosity",
                "Verbosity of process information to report [0..3].",
                ArgParseArgument::INTEGER, "INT"));

        // Set default values
        setDefaultValue(parser, "minScore",    "1000");
        setDefaultValue(parser, "nCandidates", "1");
        setDefaultValue(parser, "seedSize",    "12");
        setDefaultValue(parser, "verbosity",   "1");
            
        return parser;
    }

    void ParseArguments(int argc, char const ** argv)
    {
        result = parse(parser, argc, argv);

        // Set the ParseOk flag if ... we parsed the arguments Ok
        if (result == seqan::ArgumentParser::PARSE_OK)
        {
            parseOk = 1;
        } else {
            parseOk = 0;
        }
    }
};

SrsliParameters::SrsliParameters(int argc, char const ** argv)
{
    parser = SetupParser();
    ParseArguments(argc, argv);
    Initialize();
}
