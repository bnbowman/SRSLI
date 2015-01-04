// Author: Brett Bowman

#include <seqan/align.h>

#include "config/SeqAnConfig.hpp"
#include "utils/Align.cpp"
#include "AlignmentRecord.hpp"

namespace srsli {

    float AlignmentRecord::Accuracy()
    {
        if (accuracy < 0.1)
            accuracy = CalculateAccuracy();
        return accuracy;
    }

    // Private functions
    float AlignmentRecord::CalculateAccuracy()
    {
        return AlignmentAccuracy(Alignment, 0, 1);
    }

    // Constructors
    AlignmentRecord::AlignmentRecord(const TSegment& queryInfix,
                                     const TSegment& refInfix)
    {
        resize(rows(Alignment), 2);
        assignSource(row(Alignment, 0), queryInfix);
        assignSource(row(Alignment, 1), refInfix);
    }
}
