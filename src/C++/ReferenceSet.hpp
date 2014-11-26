// Author: Brett Bowman

#pragma once

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "SparseAlignment.hpp"

using namespace seqan;

namespace srsli {

    class ReferenceSet {

    private:
        std::string filename;
        std::string faiFilename;
        FaiIndex faiIndex;
        StringSet<CharString> ids;
        StringSet<Dna5String> seqs;
        size_t size;
        size_t seqCount;

    public:
        size_t Size() const;
        size_t Length() const;
        StringSet<CharString> Ids() const;
        StringSet<Dna5String> Sequences() const;
        seqan::FaiIndex FaiIndex() const;

    public:
        template<typename TConfig = FindSeedsConfig<>>
        Index<StringSet<Dna5String>, typename TConfig::IndexType> 
            GetIndex();

    public:
        ReferenceSet(const std::string& filename_);
    };
}
#include "ReferenceSetImpl.hpp"
