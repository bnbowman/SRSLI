// Author: David Alexander

#pragma once

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

namespace srsli {

    /// A struct for storing records as they are read
    struct SequenceRecord {
        CharString Id;
        Dna5String Seq;
        CharString Qual;
    };

    /// brief A pairwise alignment
    class SequenceReader {

    private:
        std::string Filename;
        SequenceRecord Record;
        SequenceStream Stream;
        size_t CurrentIdx;

    public:
        // Read the next record in the stream, if any
        bool GetNext(std::pair<size_t, SequenceRecord>&);

    public:
        SequenceReader(const std::string& filename);
    };
}