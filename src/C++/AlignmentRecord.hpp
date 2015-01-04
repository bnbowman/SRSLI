// Author: Brett Bowman

#pragma once

#include "config/SeqAnConfig.hpp"

namespace srsli {

    class AlignmentRecord {

    private:
        float accuracy;

    private:
        float CalculateAccuracy();

    public:
        TAlign Alignment;
        float Accuracy();
        //size_t QueryStart;
        //size_t QueryEnd;
        //size_t ReferenceStart;
        //size_t ReferenceEnd;
        //size_t Length;

    public:
        AlignmentRecord(const TSegment& queryInfix,
                        const TSegment& refInfix);
    };
}
