#include <string>
#include <vector>

#pragma once

#define API_MAJOR 0
#define API_MINOR 1
#define API_PATCH 0
#define API_DATE "November 10th, 2014"

namespace srsli
{
    class Version
    {
    public:
        static int Major();
        static int Minor();
        static int Patch();
        static std::string Date();

        // Sadly SWIG doesn't support boost::tuple
        static std::vector<int> VersionTuple();

        static std::string VersionString();
    };
}