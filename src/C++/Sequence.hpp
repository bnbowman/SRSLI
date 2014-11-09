// Author: David Alexander

#pragma once

#include <string>

namespace SRSLI
{
    char ComplementaryBase(char base);
    std::string Complement(const std::string& input);
    std::string Reverse(const std::string& input);
    std::string ReverseComplement(const std::string& input);
}