// Author: David Alexander, Brett Bowman

#pragma once

#include <string>
#include <vector>

namespace srsli {
    /// \brief A pairwise alignment
    class PairwiseAlignment {
    private:
        std::string target_;
        std::string query_;
        std::string transcript_;

    public:
        // target string, including gaps; usually the "reference"
        std::string Target() const;

        // query string, including gaps; usually the "read"
        std::string Query() const;

        // transcript as defined by Gusfield pg 215.
        std::string Transcript() const;

    public:
        float Accuracy() const;
        int Matches() const;
        int Errors() const;
        int Mismatches() const;
        int Insertions() const;
        int Deletions() const;
        int Length() const;

    public:
        PairwiseAlignment(const std::string& target,
                const std::string& query);
    };

    struct SmithWatermanParams {
        float MatchScore;
        float MismatchScore;
        float InsertScore;
        float DeleteScore;
        float BranchScore;
        float MergeScore;

        SmithWatermanParams(
                float matchScore,
                float mismatchScore,
                float insertScore,
                float deleteScore,
                float branchScore,
                float mergeScore);
    };

    SmithWatermanParams DefaultSmithWatermanParams();

    PairwiseAlignment* Align(const std::string& target,
            const std::string& query,
            SmithWatermanParams params = DefaultSmithWatermanParams()); // NOLINT

    // These calls return an array, same len as target, containing indices into the query string.
    std::vector<int> TargetToQueryPositions(const std::string& transcript);
    std::vector<int> TargetToQueryPositions(const PairwiseAlignment& aln);
}