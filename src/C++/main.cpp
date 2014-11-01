// Author: Brett Bowman

#include <boost/program_options.hpp>
#include <seqan/sequence.h>

#include "ArgumentParser.hpp"
#include "Version.hpp"


using namespace std;
using namespace seqan;
using namespace boost;
using namespace srsli;

namespace po = boost::program_options;


// Entry point
int main(int argc, char *argv[])
{
    string input;
    string reference;

    po::options_description options_desc("Allowed options");
    options_desc.add_options()
            ("query,q", po::value<string>(&input)->default_value(""), "Query sequences in FASTA format.")
            ("reference,r", po::value<string>(&input)->default_value(""), "Reference sequences in FASTA format.")
            ;

    // Exit if arguments cannot be parsed
    if (ArgumentParser::parseArguments(argc, argv, options_desc)) return 1;

    // Show help, if there is no input or input does not exist
    if (!input.compare("") || !ArgumentParser::existsTest(input))
    {
        ArgumentParser::usage(options_desc, "Input file does not exist!");
    }
}