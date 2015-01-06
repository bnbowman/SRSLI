// Author: Brett Bowman

#include <stdio.h>
#include <seqan/align.h>

#include "config/SeqAnConfig.hpp"
#include "utils/Align.cpp"
#include "utils/RegionT.cpp"
#include "AlignmentRecord.hpp"

namespace srsli {

    // Public Functions
    float AlignmentRecord::Accuracy()
    {
        EnsureAlignmentIsClipped();
        if (accuracy < 0.1)
            accuracy = AlignmentAccuracy(Alignment, 0, 1);
        return accuracy;
    }

    int AlignmentRecord::UnalignedQueryStartBases()
    {
        EnsureAlignmentIsClipped();
        if (unalignedQueryStartBases < 0)
            unalignedQueryStartBases = CountUnalignedStartBases(QueryRow());
        return unalignedQueryStartBases;
    }

    int AlignmentRecord::UnalignedQueryEndBases()
    {
        EnsureAlignmentIsClipped();
        if (unalignedQueryEndBases < 0)
            unalignedQueryEndBases = CountUnalignedEndBases(QueryRow());
        return unalignedQueryEndBases;
    }

    int AlignmentRecord::UnalignedReferenceStartBases()
    {
        EnsureAlignmentIsClipped();
        if (unalignedReferenceStartBases < 0)
            unalignedReferenceStartBases = CountUnalignedStartBases(ReferenceRow());
        return unalignedReferenceStartBases;
    }

    int AlignmentRecord::UnalignedReferenceEndBases()
    {
        EnsureAlignmentIsClipped();
        if (unalignedReferenceEndBases < 0)
            unalignedReferenceEndBases = CountUnalignedEndBases(ReferenceRow());
        return unalignedReferenceEndBases;
    }

    size_t AlignmentRecord::QueryStart()
    {
        return AlignmentRegion.queryStart + UnalignedQueryStartBases();
    }

    size_t AlignmentRecord::QueryEnd()
    {
        return AlignmentRegion.queryEnd - UnalignedQueryEndBases();
    }

    size_t AlignmentRecord::ReferenceStart()
    {
        return AlignmentRegion.refStart + UnalignedReferenceStartBases();
    }

    size_t AlignmentRecord::ReferenceEnd()
    {
        return AlignmentRegion.refEnd - UnalignedReferenceEndBases();
    }

    std::string AlignmentRecord::toM1Record()
    {
        char buffer [200];
        sprintf(buffer, "%s %s %s %s %ld %.4f %zu %zu %zu %zu %zu %zu %zu", 
                        "QueryName", "ReferenceName", "0", "0", 
                        Score, Accuracy(),
                        ReferenceStart(), ReferenceEnd(), ReferenceLength,
                        QueryStart(), QueryEnd(), QueryLength,
                        length(QueryRow()));
        return buffer;
    }

    // Private functions
    void AlignmentRecord::EnsureAlignmentIsClipped()
    {
        if (!hasBeenClipped)
        {
            ClipAlignment(Alignment);
            hasBeenClipped = true;
        }
    }

    void AlignmentRecord::Initialize()
    {
        accuracy = 0.0;
        Score = 0L;
        unalignedQueryStartBases = -1;
        unalignedQueryEndBases = -1;
        unalignedReferenceStartBases = -1;
        unalignedReferenceEndBases = -1;
        hasBeenClipped = false;
    }

    // Constructors
    AlignmentRecord::AlignmentRecord(const TDna& querySeq,
                                     const TDna& refSeq,
                                     const region_t& alignmentRegion)
    {
        // Initialize all basic variables
        Initialize();
        AlignmentRegion = alignmentRegion;
        QueryLength = length(querySeq);
        ReferenceLength = length(refSeq);

        // Infix the supplied sequences to only the required region
        TSegment queryInfix = RegionToInfix(querySeq, AlignmentRegion, 'H');
        TSegment refInfix   = RegionToInfix(refSeq, AlignmentRegion, 'V');

        // Initialize the alignment and rows using the infixes
        resize(rows(Alignment), 2);
        assignSource(QueryRow(),     queryInfix);
        assignSource(ReferenceRow(), refInfix);
    }
}
