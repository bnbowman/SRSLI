// Author: Brett Bowman

#include <boost/program_options.hpp>

#include <zlib.h>
#include <stdio.h>

#include "ArgumentParser.hpp"
#include "Utils.hpp"
#include "Version.hpp"
#include "ReferenceSet.cpp"

using namespace std;
using namespace boost;
using namespace srsli;

namespace po = boost::program_options;


// Entry point
int main(int argc, char *argv[])
{
    string query;
    string reference;

    po::options_description options_desc("Allowed options");
    options_desc.add_options()
            ("query,q", po::value<string>(&query)->default_value(""), "Query sequences in FASTA format.")
            ("reference,r", po::value<string>(&reference)->default_value(""), "Reference sequences in FASTA format.")
            ;

    // Exit if arguments cannot be parsed
    if (ArgumentParser::parseArguments(argc, argv, options_desc)) return 1;

    // Show help, if there is no input or input does not exist
    if (!query.compare("") || !ArgumentParser::existsTest(query))
    {
        ArgumentParser::usage(options_desc, "Input query file does not exist!");
    } else if (!reference.compare("") || !ArgumentParser::existsTest(reference)) {
        ArgumentParser::usage(options_desc, "Input reference file does not exist!");
    }

    // Read the reference sequences into memory
    ReferenceSet refSet = ReferenceSet(reference);


    for (unsigned i = 0; i < refSet.Length(); ++i)
        std::cout << refSet.Ids()[i] << '\t' << refSet.Sequences()[i] << '\n';
}