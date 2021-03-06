
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

using namespace seqan;

// Find seeds using the index
template<typename TConfig = FindSeedsConfig<>>
void FindSeeds(std::vector<TSeedSet>& seeds,
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
