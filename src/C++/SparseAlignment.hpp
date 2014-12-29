
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

#include "SeqAnConfig.hpp"

using namespace std;
using namespace seqan;


inline
size_t SafeSubtract(size_t size, size_t k)
{
    return size > k ? size - k : 0;
}


template <typename TSeq>
bool IsHomopolymer(const TSeq& kmer)
{
    if (length(kmer) < 2)
        return false;

    Dna b = kmer[0];

    for (size_t i = 1; i < length(kmer); i++)
    {
        if (kmer[i] != b)
            return false;
    }

    return true;
}


template<typename TConfig = FindSeedsConfig<>>
void FindSeeds(SeedSet<Simple>& seeds,
               const Dna5String& seq1, const Dna5String& seq2)
{
    Index<Dna5String, typename TConfig::IndexType> index(seq1);
    Finder<Index<Dna5String, typename TConfig::IndexType>> finder(index);

    for (size_t j = 0; j < SafeSubtract(length(seq2), TConfig::Size); j++)
    {
        Infix<Dna5String>::Type kmer = infix(seq2, j, j + TConfig::Size);

        if (IsHomopolymer(kmer))
            continue;
        
        while (find(finder, kmer))
        {
            auto pos = position(finder);
            Seed<Simple> seed(getValueI1(pos), j, TConfig::Size);

            if (!addSeed(seeds, seed, 0, Merge()))
            {
                addSeed(seeds, seed, Single());
            }
        }
    }
}

typedef std::pair<int, int> SizePosPair;

template <typename T>
bool VectorSizeCompare(vector<T> a, vector<T> b)
{
    return a.size() < b.size();
}


template<typename TConfig = FindSeedsConfig<>>
void FindSeeds(map<size_t, SeedSet<Simple>>& seeds,
               Index<StringSet<Dna5String>, typename TConfig::IndexType>& index,
               const Dna5String& seq,
               double mask = 0.999)
{
    // TODO (lhepler) : the mask refers to most common indices,
    //                  not most common kmers. You should probably fix that.
    using std::sort;

    typedef Pair<size_t, size_t> Position;

    Finder<Index<StringSet<Dna5String>, typename TConfig::IndexType>> finder(index);
    vector<vector<Position>> indexHits;
    size_t end = SafeSubtract(length(seq), TConfig::Size);

    indexHits.resize(end);


    // accumulate all the found indexHits by their index in seq
    for (size_t i = 0; i < end; i++)
    {
        Infix<const Dna5String>::Type kmer = infix(seq, i, i + TConfig::Size);

        if (IsHomopolymer(kmer))
            continue;
        
        while (find(finder, kmer))
        {
            Position pos = position(finder);
            indexHits[i].push_back(pos);
        }

        clear(finder);
    }

    //std::cout << "Found positions, sorting" << std::endl;
    
    // sort indexHits by the number of hits found
    sort(indexHits.begin(), indexHits.end(), VectorSizeCompare<Position>);

    size_t maxHits = static_cast<size_t>(mask * indexHits.size());

    // cutoff the top (1-mask) hits by index (not top (1-mask) kmers)
    for (size_t i = 0; i < maxHits; i++)
    {
        for (const Position& pos : indexHits[i])
        {
            size_t idx = getValueI1(pos);

            Seed<Simple> seed(i, getValueI2(pos), TConfig::Size);

            if (!addSeed(seeds[idx], seed, 0, Merge()))
            {
                addSeed(seeds[idx], seed, Single());
            }
        }
    }
}

// Find seeds using a qgramFinder
template<typename TConfig = FindSeedsConfig<>>
void FindSeeds2(SeedSet<Simple>& seeds,
                Index<StringSet<Dna5String>, typename TConfig::IndexType>& index,
                const Dna5String& query)
{
    typedef Shape<Dna5, typename TConfig::ShapeType> TShape;
    typedef StringSet<Dna5String> TStringSet;
    typedef Index<TStringSet, typename TConfig::IndexType> TIndex;
    typedef Finder<TIndex> TFinder;
    typedef Infix<const Dna5String> TInfix;
    typedef Iterator<const Dna5String, Standard>::Type TIterator;

    TFinder qgramFinder(index);
    for (TIterator it = begin(query, Standard()); it != end(query, Standard()) - 12; ++it)
    {
        int qPos = position(it, query);
        Dna5String qgram = infix(query, it, it+12);
                
        while(find(qgramFinder, qgram))
        {
            // Find the reference position, 
            int refPos = getValueI2( position(qgramFinder) );
            Seed<Simple> seed(qPos, refPos, TConfig::Size);
            
            if (!addSeed(seeds, seed, 0, Merge()))
            {
                addSeed(seeds, seed, Single());
            }
        }

        // Clear finder for next search.
        clear(qgramFinder); 
    }
}


// Find seeds using the index
template<typename TConfig = FindSeedsConfig<>>
void FindSeeds3(vector<TSeedSet>& seeds,
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


template<typename TConfig = FindSeedsConfig<>>
TSeedString TrimSeedChain(const TSeedString& chain,
                          unsigned short maxDiagDiff)
{
    TSeedString bestChain;
    TSeedString currentChain;
    bool prevDiagSet = false;
    int prevDiag;
    int currDiag;
    int diagDiff;

    for (auto it = begin(chain, seqan::Standard()); it != end(chain, seqan::Standard()); ++it)
    {
        currDiag = lowerDiagonal(*it);

        // If there is a previous diagonal, check the difference from current;
        if (prevDiagSet) 
        {
            diagDiff = std::abs(currDiag - prevDiag);

            // If we have a particularly large gap, stop adding to the current chain...
            if (diagDiff > maxDiagDiff) 
            { 
                bestChain = (length(currentChain) > length(bestChain)) ? currentChain : bestChain;
                clear(currentChain);
            // Otherwise append it to the current chain
            } else {
                append(currentChain, *it);
            }

        // If there was no previous diagonal, we just save the current seed
        } else {
            append(currentChain, *it);
            prevDiagSet = true;
        }
        // Whatever the location of the current/previous seeds, keep the current diagnoal
        prevDiag = currDiag;
    }

    if (length(bestChain) == 0)
        return currentChain;
    return bestChain;
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
