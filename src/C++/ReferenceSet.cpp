// Author: Brett Bowman

#include <fstream>
#include <iostream>

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "ReferenceSet.hpp"

using namespace seqan;

namespace srsli {

    size_t ReferenceSet::Size() const
    {
        return size;
    }
    size_t ReferenceSet::Length() const
    {
        return seqCount;
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

        // Build an FAI index if one doesn't exist
        //int res = build(faiIndex, filename.c_str(), faiFilename.c_str());
        //if (res != 0)
        //    throw std::runtime_error("ERROR: Could not build the index!\n");

        // Open file and create RecordReader.
        std::fstream in(filename.c_str(), std::ios::binary | std::ios::in);
        RecordReader<std::fstream, SinglePass<>> reader(in);

        // Iterate over the records in the file
        if (read2(ids, seqs, reader, seqan::Fasta()) != 0)
            throw std::runtime_error("Invalid Fasta file");

        // Set seqCount the current (non-RC'd) number of sequences
        seqCount = length(seqs);

        // Iterate over the StringSet, adding the RC sequences and sum-ing the lengths
        size = 0;
        for (size_t i = 0; i < seqCount; ++i)
        {
            Dna5String rcSeq = seqs[i];
            reverseComplement(rcSeq);
            appendValue(seqs, rcSeq);
            size += 2*length(seqs[i]);  // 2x for Forward + RC
        }
    }
}
