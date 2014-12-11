// Author: Brett Bowman

#pragma once

using namespace srsli;

class SrsliParameters {
    public:
        float maxNetIndelRate;

    void Initialize() 
    {
        maxNetIndelRate = 1.30;
    }

    SrsliParameters()
    {
        Initialize();
    }
};
