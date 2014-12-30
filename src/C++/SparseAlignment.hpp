
#pragma once

#include <iostream>
#include <stdbool.h>
#include <math.h>
#include <map>
#include <vector>
#include <typeinfo>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include "config/SeqAnConfig.hpp"
#include "config/Types.hpp"
#include "ReferenceSet.hpp"

using namespace std;
using namespace seqan;

// Find seeds using the index
template<typename TConfig = FindSeedsConfig<>>
void FindSeeds(vector<TSeedSet>& seeds,
                Index<StringSet<Dna5String>, typename TConfig::IndexType>& index,
                const size_t& refSize,
                const Dna5String& query)
{
    typedef Shape<Dna5, typename TConfig::ShapeType> TShape;
    typedef StringSet<Dna5String> TStringSet;
    typedef Index<TStringSet, typename TConfig::IndexType> TIndex;
    typedef Finder<TIndex> TFinder;
    typedef Infix<const Dna5String> TInfix;
    typedef Iterator<const Dna5String, Standard>::Type TIterator;

    TShape& shape = indexShape(index);
    hashInit(shape, begin(query, Standard()));
    for (TIterator it = begin(query, Standard()); it != end(query, Standard()) - 12; ++it)
    {
        // Hash the current Qgram, then get it's query position and number of hits
        hashNext(shape, it);
        size_t qPos  = position(it, query);
        auto hits    = getOccurrences(index, shape);
        size_t count = length(hits);

        // Skip this iteration if the Kmer doesn't exist in the reference
        if (count == 0)
            continue;

        // Compute the score for the seed based on it's frequency in the reference
        float frequency = float(count)/float(refSize);
        float score = log( 1.0/frequency );

        for (size_t i = 0; i < count; ++i)
        {
            // Find the reference position
            size_t refSeq = getValueI1( hits[i] );
            size_t refPos = getValueI2( hits[i] );
            TSeed seed = TSeed(qPos, refPos, TConfig::Size);
            setScore(seed, score);
            //std::cout << seed << " " << seqan::score(seed) << " " << score << std::endl;
            
            if (!addSeed(seeds[refSeq], seed, 0, Merge()))
            {
                addSeed(seeds[refSeq], seed, Single());
            }
        }
    }
}


template<typename TAlignConfig = GlobalAlignConfig>
Segment<const Dna5String> SeedChainToInfix(const Dna5String& seq,
                                           const String<Seed<Simple>>& chain,
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
    std::cout << "S: " << startPos << " E: " << endPos << std::endl;
    return infix(seq, 0, endPos);
}


//Stuff
template<typename TConfig = GlobalAlignConfig>
float ScoreSeedChain(const TSeedString& chain)
{
    float score = 0.0;
    for (unsigned i = 0; i < length(chain); ++i)
        score += seqan::score(chain[i]);

    return score;
}


//TODO: Why isn't the global align config working?
template<typename TAlignConfig = GlobalAlignConfig>
Align<Dna5String, ArrayGaps> SeedsToAlignment(const Dna5String& seq1, 
                                              const Dna5String& seq2,
                                              const String<TSeed>& chain,
                                              const Score<long, Simple>& scoring)
{
    auto infix1 = SeedChainToInfix(seq1, chain, 'H');
    auto infix2 = SeedChainToInfix(seq2, chain, 'V');

    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), infix1);
    assignSource(row(alignment, 1), infix2);
    AlignConfig<false, false, false, false> globalConfig;

    std::cout << "Starting alignment of sequences" << std::endl;
    long alnScore = bandedChainAlignment(alignment, chain, scoring, globalConfig);
    std::cout << "Finishing alignment of sequences" << std::endl;

    return alignment;
}
