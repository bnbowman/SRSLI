#include <string>
#include <vector>
#include <sstream>

#include "Version.hpp"

namespace srsli
{
    int Version::Major()
    {
        return API_MAJOR;
    }

    int Version::Minor()
    {
        return API_MINOR;
    }

    int Version::Patch()
    {
        return API_PATCH;
    }

    std::string Version::Date()
    {
        return API_DATE;
    }

    std::vector<int> Version::VersionTuple()
    {
        int version[3] = { API_MAJOR, API_MINOR, API_PATCH };
        return std::vector<int>(version, version + 3);
    }

    std::string Version::VersionString()
    {
        std::ostringstream VersionStream;
        VersionStream << Major() << "." << Minor() << "." << Patch();
        return VersionStream.str();
    }
}