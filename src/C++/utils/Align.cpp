#pragma once

#include <seqan/align.h>

#include "../config/SeqAnConfig.hpp"

void ClipAlignment(TAlign& alignment,
                   const int minAlignmentAnchorSize = 6)
{
    int startAnchorLength = 0;
    int endAnchorLength = 0;
    int clipStart = -1;
    int clipEnd   = -1;

    auto &query     = row(alignment, 0);
    auto &reference = row(alignment, 1);

    // Find the first string of Anchor-length matches in the alignment
    for (size_t i = 0; i < length(reference); ++i)
    {
        if (reference[i] != '-' && reference[i] == query[i]) {
            if (startAnchorLength == 0)
                clipStart = i;
            ++startAnchorLength;
        } else {
            startAnchorLength = 0;
        }
        
        if (startAnchorLength >= minAlignmentAnchorSize)
            break;
    }

    // Find the last string of Anchor-length matches in the alignment
    for (size_t i = length(reference)-1; i > 0; --i)
    {
        char refBase = reference[i];
        if (refBase != '-' && refBase == query[i]) {
            if (endAnchorLength == 0)
                clipEnd = i+1;
            ++endAnchorLength;
        } else {
            endAnchorLength = 0;
        }
        
        if (endAnchorLength >= minAlignmentAnchorSize)
            break;
    }

    // Set the identified clipping positions
    setClippedBeginPosition(query,     clipStart);
    setClippedEndPosition(query,       clipEnd);
    setClippedBeginPosition(reference, clipStart);
    setClippedEndPosition(reference,   clipEnd);
}

float AlignmentAccuracy(const Align<Dna5String, ArrayGaps>& alignment,
                        const size_t queryIdx,
                        const size_t refIdx)
{
    size_t total = 0;
    size_t matches = 0;
    size_t mismatches = 0;
    size_t insertions = 0;
    size_t deletions = 0;

    auto &query = row(alignment, queryIdx);
    auto &reference = row(alignment, refIdx);

    for (size_t i = 0; i < length(reference); ++i)
    {
        char refBase = query[i];
        char queryBase = reference[i];
        if (refBase == '-' && queryBase != '-') {
            total += 1;
            insertions += 1;
        } else if (refBase == queryBase) {
            total += 1;
            matches += 1;
        } else if (queryBase == '-') {
            total += 1;
            deletions += 1;
        } else {
            total += 1;
            mismatches += 1;
        }
    }
    return 100.0*(float)matches/(float)total;
}

size_t CountUnalignedStartBases(const TRow& row)
{
    size_t count = 0;
    for (size_t i = 0; i < length(source(row)); ++i)
    {
        auto pos = toViewPosition(row, i);
        if (pos >= 0)
            break;
        ++count;
    }
    return count;
}

size_t CountUnalignedStartBases(const TAlign& alignment,
                                const int rowNum)
{
    return CountUnalignedStartBases(row(alignment, rowNum));
}

size_t CountUnalignedEndBases(const TRow& row)
{
    auto endPos = length(row);
    size_t count = 0;
    for (size_t i = length(source(row))-1; i > 0; --i)
    {
        auto pos = toViewPosition(row, i);
        if (pos <= endPos)
            break;
        ++count;
    }
    return count;
}

size_t CountUnalignedEndBases(const TAlign& alignment,
                              const int rowNum)
{
    return CountUnalignedEndBases(row(alignment, rowNum));
}
