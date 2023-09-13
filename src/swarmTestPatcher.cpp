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
patcher::clean() {
    _inited = false;
    _robotNum = 0;
    _edgeIDs.clear();
    _edgeMax.clear();
    _tranPrefix.trans.clear();
    // _tranSuffix[i] means transition path to ith goal
    for (int i = 0; i < _tranSuffix.size(); i++) {
        _tranSuffix[i].trans.clear();
    }
    _tranSuffix.clear();
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
patcher::validState(std::vector<int> state) {
    for (int i = 0; i < _edgeIDs.size(); i++) {
        // std::cout << "Edge ID: " << _edgeIDs[i] << ", has " << state[_edgeIDs[i]] << ", exceeds " << _edgeMax[i] << std::endl;
        if ((_edgeIDs[i] < 0) || (_edgeIDs[i] >= state.size())) {
            std::cout << "Wrong Edge ID: " << _edgeIDs[i] << std::endl;
            return false;
        } else if (state[_edgeIDs[i]] > _edgeMax[i]) {
            return false;
        } else if (intInVec(i, _eidToBeRemoved) && (state[_edgeIDs[i]] > 0)) {
            return false;
        } else if (sumOverVec(state) != _robotNum) {
            std::cout << "State Sum: " << sumOverVec(state) << ", robot num: " << _robotNum << std::endl;
            return false;
        }
    }
    return true;
}

std::vector<std::pair<int, int>>
patcher::localizePatch() {
    std::vector<std::pair<int, int>> ret;
    if (_edgeIDs.size() != _edgeMax.size()) {
        std::cout << "Wrong modification size: " << _edgeIDs.size() << ", " << _edgeMax.size() << std::endl;
        return ret;
    }
    // check for prefix
    for (int i = 0; i < _tranPrefix.trans.size(); i++) {
        if (!(validState(_tranPrefix.trans[i]))) {
            ret.push_back(std::make_pair(-1, i));
        }
    }
    // and suffix
    for (int i = 0; i < _tranSuffix.size(); i++) {
        for (int j = 0; j < _tranSuffix[i].trans.size(); j++) {
            if (!(validState(_tranSuffix[i].trans[j]))) {
                ret.push_back(std::make_pair(i, j));
            }
        }
    }
    return ret;
}