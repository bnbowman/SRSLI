// Author: Brett Bowman

#include <fstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "SequenceReader.hpp"

using namespace seqan;

namespace srsli {

    // Open file and create RecordReader.
    SequenceReader::SequenceReader(const std::string& filename)
            : Filename( filename )
            , Stream( Filename.c_str() )
            , CurrentIdx( 0 )
    {
        if (!isGood(Stream))
        {
            throw std::runtime_error("ERROR: Could not open the file.");
        }
    }

    bool SequenceReader::GetNext(std::pair<size_t, SequenceRecord>& idxAndRecord)
    {
        if (atEnd(Stream))
            return false;

        size_t& idx = idxAndRecord.first;
        SequenceRecord& rec = idxAndRecord.second;

        if (readRecord(rec.Id, rec.Seq, rec.Qual, Stream) == 0)
        {
            idx = CurrentIdx++;
            return true;  // Could not record from file.
        }
        return false;
    }
}