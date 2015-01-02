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
#include "ReferenceSet.hpp"

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

float AlignmentAccuracy(const Align<Dna5String, ArrayGaps>& alignment,
                        const size_t queryIdx,
                        const size_t refIdx)
{
    size_t total = 0;
    size_t matches = 0;
    size_t mismatches = 0;
    size_t insertions = 0;
    size_t deletions = 0;

    auto query = row(alignment, queryIdx);
    auto reference = row(alignment, refIdx);

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
    //std::cout << "Match: "     << matches << std::endl;
    //std::cout << "Mismatch: "  << mismatches << std::endl;
    //std::cout << "Insertion: " << insertions << std::endl;
    //std::cout << "Deletion: "  << deletions << std::endl;
    //std::cout << "Total: "     << total << std::endl;
    return 100.0*(float)matches/(float)total;
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
                              const int maxChainBuffer)
{
    region_t alignmentRegion;
    alignmentRegion.queryStart = std::max((int)beginPositionH(chain) - maxChainBuffer, 0);
    alignmentRegion.queryEnd   = std::min((int)endPositionH(chain)   + maxChainBuffer, queryLength);
    alignmentRegion.refStart   = std::max((int)beginPositionV(chain) - maxChainBuffer, 0);
    alignmentRegion.refEnd     = std::min((int)endPositionV(chain)   + maxChainBuffer, refLength);
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

//TODO: Why isn't the global align config working?
template<typename TAlignConfig = GlobalAlignConfig>
vector<TAlign> RefChainsToAlignments(const Dna5String& querySeq, 
                                     const ReferenceSet& refSet,
                                     const vector<ReferencedSeedChain>& refChains,
                                     const Score<long, Simple>& scoring,
                                     const size_t maxAligns,
                                     const float minAccuracy,
                                     const int maxChainBuffer)
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
        auto queryInfix = RegionToInfix(querySeq, alignmentRegion, 'H');
        auto refInfix   = RegionToInfix(refSeq, alignmentRegion, 'V');

        Align<Dna5String, ArrayGaps> alignment;
        resize(rows(alignment), 2);
        std::cout << "Starting creation of infix sequences" << std::endl;
        assignSource(row(alignment, 0), queryInfix);
        std::cout << "Finished query infix sequences" << std::endl;
        assignSource(row(alignment, 1), refInfix);
        std::cout << "Finishing creation of infix sequences" << std::endl;

        std::cout << "Starting alignment of sequences" << std::endl;
        AlignConfig<false, false, false, false> globalConfig;
        long alnScore = bandedChainAlignment(alignment, shiftedString, scoring, globalConfig);
        std::cout << "Finishing alignment of sequences" << std::endl;
        alignments.push_back(alignment);

        unsigned cBeginPos = clippedBeginPosition(row(alignment, 0));
        unsigned cEndPos = clippedEndPosition(row(alignment, 0)) - 1;
        std::cout << "Aligns Seq1[" << cBeginPos << ":" << cEndPos << "]";
        std::cout << " and Seq2[" << cBeginPos << ":" << cEndPos << "]" << std::endl << std::endl;
        std::cout << length(row(alignment, 0)) << std::endl;

        float accuracy = AlignmentAccuracy(alignment, 0, 1);
        std::cout << "Accuracy: " << accuracy << std::endl;

        if (accuracy > minAccuracy) {
            std::cout << "Score = " << alnScore << std::endl;
            std::cout << alignment;
            alignments.push_back(alignment);
        } else {
            break;
        }
    }

    return alignments;
}