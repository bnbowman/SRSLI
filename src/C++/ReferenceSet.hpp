// Author: Brett Bowman

#pragma once

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "SparseAlignment.hpp"

using namespace seqan;

namespace srsli {

    class ReferenceSet {

    private:
        std::string Filename;
        StringSet<CharString> Ids_;
        StringSet<Dna5String> Seqs;

    public:
        int Length() const;
        StringSet<CharString> Ids() const;
        StringSet<Dna5String> Sequences() const;

    public:
        template<typename TConfig = FindSeedsConfig<>>
        Index<StringSet<Dna5String>, typename TConfig::IndexType> 
            GetIndex();

    public:
        ReferenceSet(const std::string& filename);
    };
}
#include "ReferenceSetImpl.hpp"
