#include <swarmTestPatcher.hpp>

using namespace swarmTest;

void
patchTrans(std::vector<std::vector<int>> &trans, std::vector<std::vector<int>> patch, int pStt, int pEnd) {
    if (pEnd < pStt) {
        std::cout << "Patch end smaller than start. Patch canceled." << std::endl;
        return;
    }
    if ((pStt < 0) || (pEnd >= trans.size())) {
        std::cout << "Patch range exceed the transition to-be-patched. Patch canceled." << std::endl;
        return;
    }
    // erase and insert is in-efficient
    trans.erase(trans.begin() + pStt, trans.begin() + pEnd + 1);
    trans.insert(trans.begin() + pStt, patch.begin(), patch.end());
}

void
patcher::makePatch(std::vector<std::vector<int>> patch, int goalID, int patchStart, int patchEnd) {
    if (goalID == -1) {
        return patchTrans(_tranPrefix.trans, patch, patchStart, patchEnd);
    } else if ((goalID >= 0) && (goalID < _tranSuffix.size())) {
        return patchTrans(_tranSuffix[goalID].trans, patch, patchStart, patchEnd);
    } else {
        std::cout << "Invalid goalID: " << goalID << std::endl;
        return;
    }
}

bool
validState(std::vector<int> state, std::vector<int> edgeIDs, std::vector<int> newMaxs) {
    for (int i = 0; i < edgeIDs.size(); i++) {
        if ((edgeIDs[i] < 0) || (edgeIDs[i] >= state.size())) {
            std::cout << "Wrong Edge ID: " << edgeIDs[i] << std::endl;
            return false;
        } else if (state[edgeIDs[i]] > newMaxs[i]) {
            return false;
        }
    }
    return true;
}

std::vector<std::pair<int, int>>
patcher::localizePatchModMax(std::vector<int> edgeIDs, std::vector<int> newMaxs) {
    std::vector<std::pair<int, int>> ret;
    if (edgeIDs.size() != newMaxs.size()) {
        std::cout << "Wrong modification size: " << edgeIDs.size() << ", " << newMaxs.size() << std::endl;
        return ret;
    }
    // check for prefix
    for (int i = 0; i < _tranPrefix.trans.size(); i++) {
        if (!(validState(_tranPrefix.trans[i], edgeIDs, newMaxs))) {
            ret.push_back(std::make_pair(-1, i));
        }
    }
    // and suffix
    for (int i = 0; i < _tranSuffix.size(); i++) {
        for (int j = 0; j < _tranSuffix[i].trans.size(); j++) {
            if (!(validState(_tranSuffix[i].trans[j], edgeIDs, newMaxs))) {
                ret.push_back(std::make_pair(i, j));
            }
        }
    }
    return ret;
}

std::vector<std::pair<int, int>>
patcher::localizePatchRmvEdg(int edgeID) {
    std::vector<int> eids, maxs;
    eids.push_back(edgeID);
    maxs.push_back(0);
    return localizePatchModMax(eids, maxs);
}