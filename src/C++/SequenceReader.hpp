// Author: David Alexander

#pragma once

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

namespace srsli {

    /// A struct for storing records as they are read
    struct SequenceRecord {
        seqan::CharString id;
        seqan::Dna5String seq;
        seqan::CharString qual;
    };

    /// brief A pairwise alignment
    class SequenceReader {

    private:
        std::string _filename;
        SequenceRecord _record;
        seqan::RecordReader<std::fstream, seqan::SinglePass<>> _reader;

    public:
        // Read the next record in the stream, if any
        int next();

    public:
        SequenceReader(const std::string& filename);
    };
}