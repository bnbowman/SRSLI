// Author: Brett Bowman

#include "ArgumentParser.hpp"

#include <iostream>
#include <sys/stat.h>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace srsli
{
    bool ArgumentParser::parseArguments(int argc, char *argv[], po::options_description &options_desc)
    {
        // Parse input options
        po::variables_map options;
        try
        {
            po::store(po::parse_command_line(argc, argv, options_desc), options);
            po::notify(options);
            return 0;
        }
        catch (std::exception &e)
        {
            std::cerr << "error: " << e.what() << "\n";
            return 1;
        }
    }

    void ArgumentParser::usage(const po::options_description &options_desc, const std::string &error)
    {
        std::cerr << "ERROR: " << error << "\n";
        std::cerr << "Usage: srsli [options]\n\n" << options_desc << "\n";
        std::cerr << "Get the latest version at https://github.com/bnbowman/srsli\n";

        exit(1);
    }

    bool ArgumentParser::existsTest(const std::string &name)
    {
        struct stat buffer;
        return (stat (name.c_str(), &buffer) == 0);
    }
}