// Author: Brett Bowman

#pragma once

#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

namespace srsli {

    template<typename TConfig>
    Index<StringSet<Dna5String>, typename TConfig::IndexType> ReferenceSet::GetIndex()
    {
        return Index<StringSet<Dna5String>, typename TConfig::IndexType>(seqs);
    }
}
