// Author: Brett Bowman

#pragma once

#include <boost/program_options.hpp>

namespace srsli
{
/**
* Helps to parse input arguments
*/
    class ArgumentParser
    {

    public:     // static constants
        /**
        * Parses command-line arguments.
        * @param[in] argc number of command-line inputs
        * @param[in] argv command-line input
        * @param[in] options_desc Boost options_description object reference
        * @return true if arguments cannot be parsed.
        */
        static bool parseArguments(
                int argc,
                char *argv[],
                boost::program_options::options_description &options_desc);
        /**
        * Outputs help for command line options
        * @param options_desc Boost options_description object reference
        * @param error The error output message
        */
        static void usage(
                const boost::program_options::options_description &options_desc,
                const std::string &error);
        /**
        * Checks if file exists.
        * @param  name Path of the file
        * @return bool
        */
        static bool existsTest(const std::string &name);

    public:     // structors
        /// Remove default constructor, static only access
        ArgumentParser() = delete;
        /// Move constructor
        ArgumentParser(ArgumentParser &&src) = delete;
        /// Copy constructor
        ArgumentParser(const ArgumentParser &src) = delete;
        /// Move assignment constructor
        ArgumentParser& operator=(ArgumentParser&& rhs) noexcept = delete;
        /// Copy assignment constructor
        ArgumentParser& operator=(const ArgumentParser& rhs) = delete;
        /// Destructor
        ~ArgumentParser() = default;
    };
}