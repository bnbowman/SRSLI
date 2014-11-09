// Author: Brett Bowman

#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "ReferenceSet.hpp"

using namespace seqan;

namespace srsli {

    int ReferenceSet::Length() const
    {
        return length(Ids_);
    }
    StringSet<CharString> ReferenceSet::Ids() const
    {
        return Ids_;
    }
    StringSet<Dna5String> ReferenceSet::Sequences() const
    {
        return Seqs;
    }

    // When initialized, read the file into memory
    ReferenceSet::ReferenceSet(const std::string& filename) : Filename( filename )
    {
        // Open file and create RecordReader.
        std::fstream in(Filename.c_str(), std::ios::binary | std::ios::in);
        RecordReader<std::fstream, SinglePass<>> Reader(in);

        // Iterate over the records in the file
        if (read2(Ids_, Seqs, Reader, seqan::Fasta()) != 0)
            throw std::runtime_error("Invalid Fasta file");
    }
}