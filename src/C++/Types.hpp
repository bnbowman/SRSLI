// Author: David Alexander

#pragma once

#include <stdint.h>
#include <string>
#include <utility>

//
// Forward declarations
//
namespace SRSLI {
    class PairwiseAlignment;
}

//
// Exception types
//

namespace srsli {
    /// \brief Abstract base class for "error"-type exceptions.  Do
    ///        not catch these.
    class ErrorBase
    {
    public:
        virtual std::string Message() const throw() = 0;
        virtual ~ErrorBase() {}
    };


    /// \brief Abstract base class for exceptions, which user code
    ///        may safely catch.
    class ExceptionBase
    {
    public:
        virtual std::string Message() const throw() = 0;
        virtual ~ExceptionBase() {}
    };


    /// \brief An exception signaling an error in ConsensusCore's internal logic.
    class InternalError : public ErrorBase
    {
    public:
        InternalError()
                : msg_("Internal error encountered!")
        {}

        explicit InternalError(const std::string& msg)
                : msg_(msg)
        {}

        std::string Message() const throw()
        {
            return msg_;
        }

    private:
        std::string msg_;
    };

    class InvalidInputError : public ErrorBase
    {
    public:
        InvalidInputError()
                : msg_("Invalid input!")
        {}

        explicit InvalidInputError(const std::string& msg)
                : msg_(msg)
        {}

        std::string Message() const throw()
        {
            return msg_;
        }

    private:
        std::string msg_;
    };

    class NotYetImplementedException : public ErrorBase
    {
    public:
        std::string Message() const throw()
        {
            return "Feature not yet implemented";
        }
    };
}