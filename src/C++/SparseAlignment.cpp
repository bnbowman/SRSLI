
#include <map>
#include <vector>

#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seeds.h>
#include <seqan/sequence.h>

#include "Exceptions.hpp"

using namespace std;
using namespace seqan;


// TODO (lhepler) : investigate default values other than 10
template <size_t TSize = 10, typename TShape = UngappedShape<TSize>, typename TIndex = IndexQGram<TShape>>
class FindSeedsConfig
{
public:
    typedef TIndex IndexType;
    static const int Size = TSize;
};

const FindSeedsConfig<> DefaultFindSeedsConfig();


#if 0
unsigned ContextCode(const DnaString& seq, size_t start, size_t k, bool& isHomopolymer)
{
    unsigned code = 0;
    unsigned b = ordValue(seq[start]);
    isHomopolymer = k > 1;

    for (size_t i = start + 1; i < start + k; i++)
    {
        unsigned c = ordValue(seq[i]);
        code += c << (2 * (i - start));
        isHomopolymer &= c == b;
    }

    return code; 
}
#endif


size_t SafeSubtract(size_t size, size_t k)
{
    return size > k ? size - k : 0;
}


bool IsHomopolymer(const Infix<DnaString>::Type& kmer)
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
               const DnaString& seq1, const DnaString& seq2,
               const TConfig& _ = FindSeedsConfig<>())
{
    Index<DnaString, typename TConfig::IndexType> index(seq1);
    Finder<Index<DnaString, typename TConfig::IndexType>> finder(index);

    for (size_t j = 0; j < SafeSubtract(length(seq2), TConfig::Size); j++)
    {
        Infix<DnaString>::Type kmer = infix(seq2, j, j + TConfig::Size);

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


template <typename T>
bool VectorSizeCompare(vector<T> a, vector<T> b)
{
    return a.size() < b.size();
}


template<typename TConfig = FindSeedsConfig<>>
void FindSeeds(map<size_t, SeedSet<Simple>>& seeds,
               const Index<DnaString, typename TConfig::IndexType>& index,
               const DnaString& seq,
               double mask = 0.999,
               const TConfig& _ = FindSeedsConfig<>())
{
    // TODO (lhepler) : the mask refers to most common indices,
    //                  not most common kmers. You should probably fix that.
    using std::sort;

    typedef Pair<size_t, size_t> Position;

    Finder<Index<DnaString, typename TConfig::IndexType>> finder(index);
    vector<vector<Position>> indexHits;
    size_t end = SafeSubtract(length(seq), TConfig::Size);

    indexHits.resize(end);

    // accumulate all the found indexHits by their index in seq
    for (size_t i = 0; i < end; i++)
    {
        Infix<DnaString>::Type kmer = infix(seq, i, i + TConfig::Size);

        if (IsHomopolymer(kmer))
            continue;
        
        while (find(finder, kmer))
        {
            Position pos = position(finder);
            indexHits[i].push_back(pos);
        }

        clear(finder);
    }

    // sort indexHits by the number of hits found
    sort(indexHits.begin(), indexHits.end(), VectorSizeCompare<Position>);

    // cutoff the top (1-mask) hits by index (not top (1-mask) kmers)
    for (size_t i = 0; i < static_cast<size_t>(mask * indexHits.size()); i++)
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

template<typename TAlignConfig>
Align<DnaString, ArrayGaps> SeedsToAlignment(const DnaString& seq1, const DnaString& seq2,
                                             const SeedSet<Simple>& seeds,
                                             const Score<int, Simple>& scoring,
                                             const TAlignConfig& config)
{
    String<Seed<Simple>> chain;
    chainSeedsGlobally(chain, seeds, SparseChaining());

    Align<DnaString, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), seq1);
    assignSource(row(alignment, 1), seq2);

    bandedChainAlignment(alignment, chain, scoring, config);

    return alignment;
}
