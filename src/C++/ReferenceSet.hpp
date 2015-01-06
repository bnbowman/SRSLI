// Author: Brett Bowman

#pragma once

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include "config/SeqAnConfig.hpp"
#include "config/Types.hpp"
#include "SparseAlignment.hpp"

using namespace seqan;

namespace srsli {

    class ReferenceSet {

    private:
        std::string filename;
        std::string faiFilename;
        FaiIndex faiIndex;
        StringSet<CharString> ids;
        StringSet<TDna> seqs;
        size_t size;
        size_t seqCount;

    public:
        vector<ReferenceRecord> Records;
        size_t Size() const;
        size_t Length() const;
        StringSet<CharString> Ids() const;
        StringSet<TDna> Sequences() const;
        seqan::FaiIndex FaiIndex() const;

    public:
        template<typename TConfig = FindSeedsConfig<>>
        Index<StringSet<TDna>, typename TConfig::IndexType> 
            GetIndex();

    public:
        ReferenceSet(const std::string& filename_);
    };
}
#include "ReferenceSetImpl.hpp"
