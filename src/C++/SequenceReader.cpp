// Author: Brett Bowman

#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "SequenceReader.hpp"

using namespace seqan;

namespace srsli {

    // Open file and create RecordReader.
    SequenceReader::SequenceReader(const std::string& filename) : _filename( filename )
    {
        std::fstream in(_filename, std::ios::binary | std::ios::in);
        seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);
        auto reader = new seqan::RecordReader<std::fstream, seqan::SinglePass<> >(in);
    }

    //int SequenceReader::next()
    //{
        // Read file record-wise.
        seqan::CharString id;
        seqan::Dna5String seq;
        seqan::CharString qual;

        while (!atEnd(reader))
        {
            if (readRecord(id, seq, qual, _reader, seqan::Fastq()) != 0)
                return 1;  // Could not record from file.

            std::cout << id << "\t" << seq << "\n";
        }
        return 0;
    }
}