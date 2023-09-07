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
      void clean();

      void addTransitionPrefix(std::vector<int> tran) {
          if (!_inited) {
              _inited = true;
              _robotNum = sumOverVec(tran);
          }
          _tranPrefix.trans.push_back(tran);
          int tranSum = sumOverVec(tran);
          if (tranSum != _robotNum) {
              std::cout << "Transition Added to Prefix has different robot number." << std::endl;
          }
      }

        void newGoalSuffix(){
            t2g newEmpty;
            _tranSuffix.push_back(newEmpty);
        }

        // add or remove a certain edge in the existing patch
        // add: add an empty assignment at each transition after the edgeID
        // remove: remove the assignment at each transition at the edgeID
        void updateEdge(bool remove, int edgeID) {
            if (!_inited) {
                std::cout << "No any transitions in the patcher yet..." << std::endl;
            }
            // modify the prefix
            for (int i = 0; i < _tranPrefix.trans.size(); i++) {
                updateEdgeOnATran(remove, edgeID, _tranPrefix.trans[i]);
            }
            // modify the suffix
            for (int j = 0; j < _tranSuffix.size(); j++) {
                for (int i = 0; i < _tranSuffix[j].trans.size(); i++) {
                    updateEdgeOnATran(remove, edgeID, _tranSuffix[j].trans[i]);
                }
            }
        }

        void updateCapacity(std::vector<int> edgeIDs, std::vector<int> edgeMax) {
            _edgeIDs = edgeIDs;
            _edgeMax = edgeMax;
        }

        // if goalID == -1, then the transition is added to the prefix part
        void addTransitionSuffix(std::vector<int> tran, int goalID){
            if (!_inited) {
                _inited = true;
                _robotNum = sumOverVec(tran);
            }
            if (goalID == -1){
                return addTransitionPrefix(tran);
            } else if ((goalID >= 0) && (goalID < _tranSuffix.size())){
                _tranSuffix[goalID].trans.push_back(tran);
                int tranSum = sumOverVec(tran);
                if (tranSum != _robotNum) {
                    std::cout << "Transition Added to Prefix has different robot number." << std::endl;
                }
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

        std::vector<std::vector<int>> getAllTransitions() {
            std::vector<std::vector<int>> ret;
            int nP = _tranPrefix.trans.size();
            for (int i = 0; i < nP; i++) {
                ret.push_back(_tranPrefix.trans[i]);
            }
            int nG = _tranSuffix.size();
            for (int g = 0; g < nG; g++) {
                int nS = _tranSuffix[g].trans.size();
                for (int i = 0; i < nS; i++) {
                    ret.push_back(_tranSuffix[g].trans[i]);
                }
            }
            return ret;
        }

        int getPrefixSize() {
            if (!_inited) {
                std::cout << "Patcher is not inited. Return -1 as prefix size. " << std::endl;
                return -1;
            }
            return _tranPrefix.trans.size();
        }

        int getSuffixSize(int goalNum) {
            if (!_inited) {
                std::cout << "Patcher is not inited. Return -1 as suffix size. " << std::endl;
                return -1;
            }
            if (goalNum == -1) {
                return getPrefixSize();
            }
            return _tranSuffix[goalNum].trans.size();
        }

        // replace _tranSuffix[goalID][patchStart,..,patchEnd] with patch.  Similarly, goalID == -1 makes it patch the prefix
        void makePatch(std::vector<std::vector<int>> patch, int goalID, int patchStart, int patchEnd);

        // return a series of locations in the original patch that require patching.  The location is a [goalID, tranID] pair
        // modify the max size of a series of edges in the map
        // std::vector<std::pair<int, int>> localizePatchModMax(std::vector<int> edgeIDs, std::vector<int> newMaxs);

        // remove a certain edge
        // std::vector<std::pair<int, int>> localizePatchRmvEdg(int edgeID);

        // general patch localizer that check
        // 1.) if the subswarm size exceeds edge max capacity.
        // 2.) if the transition has different total capacity than "robotNum"
        std::vector<std::pair<int, int>> localizePatch();

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
      bool _inited = false;

      int _robotNum = 0;

      std::vector<int> _edgeIDs, _edgeMax;

      t2g _tranPrefix;

      // _tranSuffix[i] means transition path to ith goal
      std::vector<t2g> _tranSuffix;

      void updateEdgeOnATran(bool remove, int edgeID, std::vector<int> &tran) {
          if (remove) {
              if ((edgeID >= tran.size())) {
                  std::cout << "Edge to be removed not existed" << std::endl;
                  return;
              }
              tran.erase(tran.begin() + edgeID, tran.begin() + edgeID + 1);
          } else {
              if ((edgeID > tran.size())) {
                  std::cout << "Edge to be added not existed" << std::endl;
                  return;
              }
              tran.insert(tran.begin() + edgeID + 1, 0);
          }
      }

      int
      sumOverVec(std::vector<int> v) {
          int sum = 0;
          for (int i = 0; i < v.size(); i++) {
              sum += v[i];
          }
          return sum;
      }

      bool validState(std::vector<int> state);
    };
}