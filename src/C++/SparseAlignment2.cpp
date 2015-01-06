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
#include "utils/RegionT.cpp"
#include "ReferenceSet.hpp"
#include "AlignmentRecord.cpp"

using namespace seqan;
using namespace srsli;


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

//TODO: Why isn't the global align config working?
template<typename TAlignConfig = GlobalAlignConfig>
int RefChainsToAlignments(vector<AlignmentRecord>& results,
                          const Dna5String& querySeq, 
                          const ReferenceSet& refSet,
                          const vector<ReferencedSeedChain>& refChains,
                          const Score<long, Simple>& scoring,
                          const size_t maxAligns,
                          const float minAccuracy,
                          const int maxChainBuffer,
                          const int alignmentAnchorSize)
{
    AlignConfig<false, false, true, true> globalConfig;
    ReferencedSeedChain refChain;
    ReferenceRecord refRec;
    TDna* refSeqPtr;
    TSeedChain* seedChain;

    for (size_t i = 0; i < maxAligns; ++i)
    {
        refChain = refChains[i];
        size_t refIdx = refChain.referenceIndex;
        seedChain = &refChain.chain;
        refRec = refSet.Records[refIdx];
        refSeqPtr = refRec.seq;

        String<TSeed> seedString;
        for (size_t j = 0; j < seedChain->size(); ++j)
        {
            appendValue(seedString, seedChain->at(j));
        }

        region_t alignmentRegion = ChoseAlignmentRegion(seedString, 
                                                        length(querySeq), 
                                                        length(*refSeqPtr), 
                                                        maxChainBuffer);
        String<TSeed> shiftedString = ShiftSeedString(seedString, alignmentRegion);

        // Create an AlignmentRecord from the sequences and the selected region
        AlignmentRecord alnRec(querySeq, *refSeqPtr, alignmentRegion);

        std::cout << "Starting alignment of sequences" << std::endl;
        alnRec.Score = bandedChainAlignment(alnRec.Alignment, shiftedString, scoring, globalConfig);
        std::cout << "Finishing alignment of sequences" << std::endl;
        std::cout << "Accuracy: " << alnRec.Accuracy() << std::endl;

        if (alnRec.Accuracy() > minAccuracy) {
            results.push_back(alnRec);
        } else {
            break;
        }
    }
    return 0;
}
