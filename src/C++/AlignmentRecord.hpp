// Author: Brett Bowman

#pragma once

#include "config/SeqAnConfig.hpp"

namespace srsli {

    class AlignmentRecord {

    private:
        // Value variables
        float accuracy;
        int unalignedQueryStartBases;
        int unalignedQueryEndBases;
        int unalignedReferenceStartBases;
        int unalignedReferenceEndBases;

        // Flag variables
        bool hasBeenClipped;

    private:
        void Initialize();
        void EnsureAlignmentIsClipped();
        int UnalignedQueryStartBases();
        int UnalignedQueryEndBases();
        int UnalignedReferenceStartBases();
        int UnalignedReferenceEndBases();

    public:
        TAlign Alignment;
        region_t AlignmentRegion;
        size_t QueryLength;
        size_t ReferenceLength;
        long Score;

    public:
        float Accuracy();
        inline TRow& QueryRow() { return row(Alignment, 0); }
        inline TRow& ReferenceRow() { return row(Alignment, 1); }
        size_t QueryStart();
        size_t QueryEnd();
        size_t ReferenceStart();
        size_t ReferenceEnd();
        std::string toM1Record();

    public:
        AlignmentRecord(const TDna& querySeq,
                        const TDna& refSeq,
                        const region_t& alnRegion);
    };
}
