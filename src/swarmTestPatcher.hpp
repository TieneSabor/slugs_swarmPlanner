#pragma once

#include <utility>
#include <iostream>
#include <vector>

namespace swarmTest{
    typedef struct trans2Goal{
        std::vector<std::vector<int>> trans; // trans[i] is the ith step in the path to the certain goal
    }t2g;

    class patcher{
    public:
        void addTransitionPrefix(std::vector<int> tran){
            _tranPrefix.trans.push_back(tran);
        }

        void newGoalSuffix(){
            t2g newEmpty;
            _tranSuffix.push_back(newEmpty);
        }

        // if goalID == -1, then the transition is added to the prefix part
        void addTransitionSuffix(std::vector<int> tran, int goalID){
            if (goalID == -1){
                return addTransitionPrefix(tran);
            } else if ((goalID >= 0) && (goalID < _tranSuffix.size())){
                _tranSuffix[goalID].trans.push_back(tran);
            } else{
                std::cout << "Invalid goalID: " << goalID << std::endl;
                return;
            }
        }

        // if goalID == -1, then the prefix transition is returned
        std::vector<int> getTransition(int goalID, int tranID){
            if((goalID == -1)&&(tranID < _tranPrefix.trans.size())){
                return _tranPrefix.trans[tranID];
            } else if ((goalID >= 0) && (goalID < _tranSuffix.size())){
                if ((tranID >= 0)&&(tranID < _tranSuffix[goalID].trans.size())){
                    return _tranSuffix[goalID].trans[tranID];
                } else {
                    std::cout << "Invalid tranID: " << tranID << " at goal: " << goalID << std::endl;
                    std::vector<int> ret;
                    return ret;
                }
            } else{
                std::cout << "Invalid goalID: " << goalID << std::endl;
                std::vector<int> ret;
                return ret;
            }
        }

        // replace _tranSuffix[goalID][patchStart,..,patchEnd] with patch.  Similarly, goalID == -1 makes it patch the prefix
        void makePatch(std::vector<std::vector<int>> patch, int goalID, int patchStart, int patchEnd);

        // return a series of locations in the original patch that require patching.  The location is a [goalID, tranID] pair
        // modify the max size of a series of edges in the map
        std::vector<std::pair<int, int>> localizePatchModMax(std::vector<int> edgeIDs, std::vector<int> newMaxs);

        // remove a certain edge
        std::vector<std::pair<int, int>> localizePatchRmvEdg(int edgeID);

        // utilities
        void printState(std::vector<int> state){
            std::cout << "{ ";
            for (int i=0; i< state.size(); i++){
                std::cout << state[i] << ", ";
            }
            std::cout << "} ";
        }

        void printTransitions(){
            // print the prefix
            std::cout << "Prefix Transitions" << std::endl;
            for (int i=0; i<_tranPrefix.trans.size(); i++){
                std::cout << "Step [" << i << ", " << i+1 << "]: ";
                printState(_tranPrefix.trans[i]);
                std::cout << std::endl;
            }
            // suffix
            std::cout << "Suffix Transitions" << std::endl;
            for (int i=0; i<_tranSuffix.size(); i++){
                std::cout << "Goal ID: " << i << std::endl;
                for (int j=0; j<_tranSuffix[i].trans.size(); j++){
                    std::cout << "Step [" << j << ", " << j+1 << "]: ";
                    printState(_tranSuffix[i].trans[j]);
                    std::cout << std::endl;
                }
            }
        }
    private:
        t2g _tranPrefix;

        // _tranSuffix[i] means transition path to ith goal
        std::vector<t2g> _tranSuffix;
    };
}