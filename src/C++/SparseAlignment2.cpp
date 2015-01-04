#pragma once

#include <vector>
#include <algorithm>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include "config/SeqAnConfig.hpp"
#include "config/Types.hpp"
#include "utils/SeedString.cpp"
#include "utils/Align.cpp"
#include "ReferenceSet.hpp"
#include "AlignmentRecord.cpp"

using namespace seqan;
using namespace srsli;


Segment<const Dna5String> SeedChainToInfix(const Dna5String& seq,
                                           const vector<TSeed>& chain,
                                           char axis)
{
    int startPos, endPos;
    if (axis == 'H')
    {
       startPos = beginPositionH(front(chain));
       endPos   = endPositionH(back(chain));
    } else {
       startPos = beginPositionV(front(chain));
       endPos   = endPositionV(back(chain));
    }
    return infix(seq, 0, endPos);
}

struct region_t {
    int queryStart;
    int queryEnd;
    int refStart;
    int refEnd;
};

// Based on the supplied chain and buffer size, decide the regions of the sequences to align
region_t ChoseAlignmentRegion(const String<TSeed>& chain,
                              const int queryLength,
                              const int refLength,
                              const int maxChainBuffer,
                              const float netMaxIndels = 1.3)
{
    region_t alignmentRegion;
    int querySeedStart = beginPositionH(chain);
    int querySeedEnd   = endPositionH(chain);

    // First find the Query Start and end positions, given some maximum buffer size
    alignmentRegion.queryStart = std::max(querySeedStart - maxChainBuffer, 0);
    alignmentRegion.queryEnd   = std::min(querySeedEnd   + maxChainBuffer, queryLength);
    
    // Second, decide based on the Query buffer used how much sequence to take from the Reference
    int queryStartBuffer = (querySeedStart - alignmentRegion.queryStart) * netMaxIndels;
    int queryEndBuffer   = (alignmentRegion.queryEnd - querySeedEnd) * netMaxIndels;
    
    // Perform the same calculation as the for the query for the reference,
    //    using the new buffer sizes
    alignmentRegion.refStart   = std::max((int)beginPositionV(chain) - queryStartBuffer, 0);
    alignmentRegion.refEnd     = std::min((int)endPositionV(chain)   + queryEndBuffer, refLength);

    return alignmentRegion;
}

String<TSeed> ShiftSeedString(const String<TSeed>& string,
                              const region_t alignmentRegion)
{
    String<TSeed> output;
    for (size_t i = 0; i < length(string); ++i)
    {   
        TSeed newSeed(beginPositionH(string[i]) - alignmentRegion.queryStart,
                      beginPositionV(string[i]) - alignmentRegion.refStart,
                      endPositionH(string[i])   - alignmentRegion.queryStart,
                      endPositionV(string[i])   - alignmentRegion.refStart);
        appendValue(output, newSeed);
    }
    return output;
}

Segment<const Dna5String> RegionToInfix(const Dna5String& seq,
                                        region_t alignmentRegion,
                                        char axis)
{
    int startPos, endPos;
    if (axis == 'H')
    {
       startPos = alignmentRegion.queryStart;
       endPos   = alignmentRegion.queryEnd;
    } else {
       startPos = alignmentRegion.refStart;
       endPos   = alignmentRegion.refEnd;
    }
    return infix(seq, startPos, endPos);
}

void ClipAlignment(TAlign& alignment,
                   const int minAlignmentAnchorSize)
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

//TODO: Why isn't the global align config working?
template<typename TAlignConfig = GlobalAlignConfig>
vector<TAlign> RefChainsToAlignments(const Dna5String& querySeq, 
                                     const ReferenceSet& refSet,
                                     const vector<ReferencedSeedChain>& refChains,
                                     const Score<long, Simple>& scoring,
                                     const size_t maxAligns,
                                     const float minAccuracy,
                                     const int maxChainBuffer,
                                     const int alignmentAnchorSize)
{
    vector<TAlign> alignments;
    ReferencedSeedChain refChain;
    Dna5String refSeq;
    TSeedChain* seedChain;

    for (size_t i = 0; i < maxAligns; ++i)
    {
        refChain = refChains[i];
        size_t refIdx = std::get<0>(refChain);
        seedChain = &std::get<1>(refChain);
        refSeq = refSet.Sequences()[refIdx];

        String<TSeed> seedString;
        for (size_t j = 0; j < seedChain->size(); ++j)
        {
            appendValue(seedString, seedChain->at(j));
        }

        region_t alignmentRegion = ChoseAlignmentRegion(seedString, 
                                                        length(querySeq), 
                                                        length(refSeq), 
                                                        maxChainBuffer);
        String<TSeed> shiftedString = ShiftSeedString(seedString, alignmentRegion);
        TSegment queryInfix = RegionToInfix(querySeq, alignmentRegion, 'H');
        TSegment refInfix   = RegionToInfix(refSeq, alignmentRegion, 'V');

        // Create an AlignmentRecord from the sequence infixes
        AlignmentRecord alnRec(queryInfix, refInfix);

        std::cout << "Starting alignment of sequences" << std::endl;
        AlignConfig<false, false, true, true> globalConfig;
        long alnScore = bandedChainAlignment(alnRec.Alignment, shiftedString, scoring, globalConfig);
        std::cout << "Finishing alignment of sequences" << std::endl;

        std::cout << "Clipping alignment to the core matching region" << std::endl;
        ClipAlignment(alnRec.Alignment, alignmentAnchorSize);
        std::cout << alnRec.Alignment << std::endl;
        std::cout << "Finishing clipping alignment" << std::endl;
        std::cout << alnRec.Accuracy() << std::endl;

        float accuracy = AlignmentAccuracy(alnRec.Alignment, 0, 1);
        std::cout << "Accuracy: " << accuracy << std::endl;

        if (accuracy > minAccuracy) {
            //auto& row0 = row(alignment, 0);
            //auto& row1 = row(alignment, 1);
            //std::cout << CountUnalignedStartBases(row0) << std::endl;
            //std::cout << CountUnalignedStartBases(row1) << std::endl;
            //std::cout << CountUnalignedEndBases(row0) << std::endl;
            //std::cout << CountUnalignedEndBases(row1) << std::endl;
            //alignments.push_back(alignment);
        } else {
            break;
        }
    }

    return alignments;
}
