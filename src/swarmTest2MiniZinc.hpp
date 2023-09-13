#pragma once

// for c forking stuffs
#include <errno.h>
// #include <ext/stdio_filebuf.h>
#include <cstdio>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
// others
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace swarmTest {

class regionMap {
  public:
    regionMap(unsigned int regionNumber);
    // ~regionMap();

    void addEdge(int fromIndex, int toIndex);

    void rmvEdge(int fromIndex, int toIndex);

    // return a list of in-neighbors
    std::vector<int> getInRegion(int regionIndex);

    // return a list of out-neighbors
    std::vector<int> getOutRegion(int regionIndex);

    // print
    void printRegionAM();

  private:
    // Adjancy Matrix of the region map, regionAM[i][j]==1 iff edge "i to j" exists
    std::vector<std::vector<bool>> _regionAM;
};

enum printDest {
    toTerminal,
    toError,
    toFile,
    toFork,
    toPatcher
};

enum printFor {
    patching,
    reassignment,
    planassignment
};

struct clause {
    std::vector<int> beforeRegionIndex;
    std::vector<bool> beforeComplemented;
    std::vector<int> afterRegionIndex;
    std::vector<bool> afterComplemented;
};

class convert2MiniZinc {
  public:
    convert2MiniZinc(std::vector<std::string> regionNames);
    // ~convert2MiniZinc();

    int getRegionNum() {
        return _maxNs.size();
    }

    std::vector<int> getRegionCap() {
        return _maxNs;
    }

    std::vector<int> getRegionInit() {
        return _iniNs;
    }

    std::string getRegionName(int rid) {
        if ((rid >= _regionNames.size()) || (rid < 0)) {
            std::cout << "Region ID " << rid << " is not existed" << std::endl;
            return "";
        }
        return _regionNames[rid];
    }

    // if the edge does not exist, return -1
    int getEdgeIDByIndex(int fromIndex, int toIndex);

    // if the edge does not exist, return -1
    int getEdgeIDByName(std::string nameFrom, std::string nameTo);

    std::pair<int, int> regionFromToByEdgeIndex(int eid);

    // if ***N is -1, then the value remains unchanged
    void setRegionConstByIndex(int index, int maxN, int iniN, int fnlN) {
        if (maxN != -1)
            _maxNs[index] = maxN;
        if (iniN != -1)
            _iniNs[index] = iniN;
        if (fnlN != -1)
            _fnlNs[index] = fnlN;
    }

    void setRegionConstByName(std::string Name, int maxN, int iniN, int fnlN) {
        setRegionConstByIndex(_regionIDsByName[Name], maxN, iniN, fnlN);
    }

    // if directed is true, add only one edge: from -> to. Else, add both edges: from <-> to
    // if remove, remove the edge
    void addEdgeByName(bool directed, bool remove, std::string fromName, std::string toName);

    // if directed is true, add only one edge: from -> to. Else, add both edges: from <-> to
    // if remove, remove the edge
    void addEdgeByIndex(bool directed, bool remove, int fromIndex, int toIndex);

    // "A" : region A_i > 0
    // "-A": region A_i == 0
    // "A'": region A_(i+1) > 0
    void addLiteral2LastClauseByName(std::string name) {
        bool com = false;
        bool bef = true;
        if (name[0] == '-') {
            com = true;
            name = name.substr(1);
        }
        if (name.back() == '\'') {
            bef = false;
            name = name.substr(0, name.size() - 1);
        }
        int idL = _regionIDsByName[name];
        addLiteral2LastClauseByIndex(com, bef, idL);
    }

    // complemented == true -> complemented literal
    // befAft == true       -> before predicate; else, after predicate
    // idL: index of the literal
    void addLiteral2LastClauseByIndex(bool complemented, bool before, int idL) {
        addLiteral2CertainClauseByIndex(complemented, before, idL, _clauses.size() - 1);
    }

    void addLiteral2CertainClauseByIndex(bool complemented, bool before, int idL, int idC) {
        if (before) {
            _clauses[idC].beforeComplemented.push_back(complemented);
            _clauses[idC].beforeRegionIndex.push_back(idL);
        } else {
            _clauses[idC].afterComplemented.push_back(complemented);
            _clauses[idC].afterRegionIndex.push_back(idL);
        }
    }

    void newClause() {
        struct clause c;
        _clauses.push_back(c);
    }

    void clearClauses() {
        _clauses.clear();
    }

    void printMiniZinc(printDest dest, printFor reason, std::string fileName, unsigned int layerNumber, std::vector<std::vector<int>> assignment, bool finalBC);

    void printPatch2Dot(std::string fileName);

    void printPatch2Dot(std::string fileName, std::vector<std::vector<int>> patch);

    // Utilities
    void printTransition(std::vector<int> transition);

    void printClauses() {
        std::cout << "Print Stored Clauses in CNF: " << std::endl;
        int nC = _clauses.size();
        for (int i = 0; i < nC; i++) {
            int iB = _clauses[i].beforeRegionIndex.size();
            for (int j = 0; j < iB; j++) {
                std::cout << _clauses[i].beforeRegionIndex[j] << "_pre \\/ ";
            }
            int iA = _clauses[i].afterRegionIndex.size();
            for (int j = 0; j < iA; j++) {
                std::cout << _clauses[i].afterRegionIndex[j] << "_post \\/ ";
            }
            std::cout << std::endl;
        }
    }

    std::vector<int> getRegionInNeighborByIndex(int id) {
        return _rM.getInRegion(id);
    }

    std::vector<int> getRegionOutNeighborByIndex(int id) {
        return _rM.getOutRegion(id);
    }

    std::vector<std::vector<int>> getPatch() {
        return _patch;
    }

    std::vector<int> edge2RegionState(bool before, std::vector<int> edgeState) {
        int nR = _regionNames.size();
        std::vector<int> retBf(nR, 0), retAf(nR, 0);
        int edgeID = 0;
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                retBf[j] += edgeState[edgeID];             // before: sum of out edge
                retAf[outEdges[k]] += edgeState[edgeID];   // after: sum of in edge
                edgeID++;
            }
        }
        // if before, return the region state before the edge state
        if (before) {
            return retBf;
        } else {
            return retAf;
        }
    }

    std::vector<int> getEdgeMax() {
        std::vector<int> eM;
        int nR = _regionNames.size();
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                eM.push_back(std::min(_maxNs[j], _maxNs[outEdges[k]]));
            }
        }
        return eM;
    }

    void setEdgeLiterals(std::vector<std::vector<std::pair<int, bool>>> literals) {
        _edgeLiterals = literals;
    }

    void printEdgeLiterals();

    void testPatchWithEdgeLiterals();

  private:
    std::vector<std::string> _regionNames;

    std::map<std::string, int> _regionIDsByName;

    std::vector<int> _maxNs;

    std::vector<int> _iniNs;

    std::vector<int> _fnlNs;

    regionMap _rM;

    std::vector<struct clause> _clauses;

    std::vector<std::vector<std::pair<int, bool>>> _edgeLiterals;

    bool _edgeFTInited = false, _fTEdgeInited = false;
    std::vector<std::pair<int, int>> _edgeFromTo;
    std::vector<std::vector<int>> _fromToEdge;

    std::string miniZincString(unsigned int layernumber, bool finalBC);

    std::string miniZincAssign(std::vector<std::vector<int>> assignment, bool finalBC);

    std::string miniZincReassign(std::vector<std::vector<int>> assignment, bool finalBC);

    // for working with the patcher
    std::vector<std::vector<int>> _patch;

    // for working with the dot display
    std::string regionState2NodeString(std::vector<int> rS, int stateID);
};
}   // namespace swarmTest