// Author: David Alexander

#include <string>

namespace SRSLI {

    char ComplementaryBase(char base) {
        return ComplementArray[(int) base];
    }

    std::string Complement(const std::string &input) {
        std::string output(input.length(), 127);
        for (unsigned int i = 0; i < input.length(); i++) {
            output[i] = ComplementArray[(int) input[i]];
        }
        return output;
    }

    std::string Reverse(const std::string &input) {
        return std::string(input.rbegin(), input.rend());
    }

    std::string ReverseComplement(const std::string &input) {
        return Reverse(Complement(input));
    }
}