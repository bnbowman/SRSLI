// Author: Brett Bowman

#include <fstream>
#include <iostream>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "ReferenceSet.hpp"

using namespace seqan;

namespace srsli {

    int ReferenceSet::Length() const
    {
        return length(ids);
    }
    StringSet<CharString> ReferenceSet::Ids() const
    {
        return ids;
    }
    StringSet<Dna5String> ReferenceSet::Sequences() const
    {
        return seqs;
    }
    FaiIndex ReferenceSet::FaiIndex() const
    {
        return faiIndex;
    }

    // When initialized, read the file into memory
    ReferenceSet::ReferenceSet(const std::string& filename_)
            : filename( filename_ )
            , faiFilename( filename + ".fai" )
    {
        std::cout << "file: " << filename    << std::endl;
        std::cout << "fai: "  << faiFilename << std::endl;

        int res = build(faiIndex, filename.c_str(), faiFilename.c_str());
        if (res != 0)
            throw std::runtime_error("ERROR: Could not build the index!\n");

        // Open file and create RecordReader.
        std::fstream in(filename.c_str(), std::ios::binary | std::ios::in);
        RecordReader<std::fstream, SinglePass<>> reader(in);

        // Iterate over the records in the file
        if (read2(ids, seqs, reader, seqan::Fasta()) != 0)
            throw std::runtime_error("Invalid Fasta file");
    }
}
