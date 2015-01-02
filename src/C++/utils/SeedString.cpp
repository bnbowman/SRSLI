#pragma once
#include <vector>
#include <seqan/seeds.h>

#include "../config/SeqAnConfig.hpp"


template<typename TSeed>
int beginPositionV(const String<TSeed>& string)
{
    return beginPositionV( front(string) );
}

template<typename TSeed>
int endPositionV(const String<TSeed>& string)
{
    return endPositionV( back(string) );
}

template<typename TSeed>
int beginPositionH(const String<TSeed>& string)
{
    return beginPositionH( front(string) );
}

template<typename TSeed>
int endPositionH(const String<TSeed>& string)
{
    return endPositionH( back(string) );
}
// End Utility Functions
