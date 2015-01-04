#pragma once

#include <seqan/align.h>

#include "../config/SeqAnConfig.hpp"

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
