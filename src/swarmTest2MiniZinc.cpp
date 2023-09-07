#include "swarmTest2MiniZinc.hpp"

using namespace swarmTest;

regionMap::regionMap(unsigned int regionNumber) : _regionAM(regionNumber, std::vector<bool>(regionNumber, false)) {}

void
regionMap::addEdge(int fromIndex, int toIndex) {
    _regionAM[fromIndex][toIndex] = true;
}

void
regionMap::rmvEdge(int fromIndex, int toIndex) {
    _regionAM[fromIndex][toIndex] = false;
}

std::vector<int>
regionMap::getInRegion(int regionIndex) {
    std::vector<int> ret;
    int nR = _regionAM.size();
    for (int i = 0; i < nR; i++) {
        if (_regionAM[i][regionIndex] == true)
            ret.push_back(i);
    }
    return ret;
}

// return a list of out-neighbors
std::vector<int>
regionMap::getOutRegion(int regionIndex) {
    std::vector<int> ret;
    int nR = _regionAM.size();
    for (int i = 0; i < nR; i++) {
        if (_regionAM[regionIndex][i] == true)
            ret.push_back(i);
    }
    return ret;
}

// print
void
regionMap::printRegionAM() {
    int nR = _regionAM.size();
    for (int i = 0; i < nR; i++) {
        for (int j = 0; j < nR; j++) {
            if (_regionAM[i][j] == true)
                std::cout << " " << 1;
            else
                std::cout << " " << 0;
        }
        std::cout << std::endl;
    }
}

// utilities for fork stuff
int
waitPidHandy(int pid) {
    int status;
    // the parent process calls waitpid() on the child
    // waitpid() system call suspends execution of
    // calling process until a child specified by pid
    // argument has changed state
    // see wait() man page for all the flags or options
    // used here
    if (waitpid(pid, &status, 0) > 0) {

        if (WIFEXITED(status) && !WEXITSTATUS(status)) {
            printf("program execution successful\n");
            return 0;
        }

        else if (WIFEXITED(status) && WEXITSTATUS(status)) {
            if (WEXITSTATUS(status) == 127) {
                // execv failed
                printf("execv failed\n");
                return -1;
            } else {
                printf("program terminated normally,"
                       " but returned a non-zero status\n");
                return status;
            }
        } else {
            printf("program didn't terminate normally\n");
            return -1;
        }
    } else {
        // waitpid() failed
        printf("waitpid() failed\n");
        return -1;
    }
}

// convert2MiniZinc implementations
convert2MiniZinc::convert2MiniZinc(std::vector<std::string> regionNames)
    : _rM(regionNames.size()), _maxNs(regionNames.size()), _iniNs(regionNames.size()), _fnlNs(regionNames.size()) {
    // create lists for name-index pairs
    int nR = regionNames.size();
    for (int i = 0; i < nR; i++) {
        _regionIDsByName[regionNames[i]] = i;
        _regionNames.push_back(regionNames[i]);
    }
}

void
convert2MiniZinc::addEdgeByName(bool directed, bool remove, std::string fromName, std::string toName) {
    int fromIndex = _regionIDsByName[fromName];
    int toIndex = _regionIDsByName[toName];
    addEdgeByIndex(directed, remove, fromIndex, toIndex);
}

void
convert2MiniZinc::addEdgeByIndex(bool directed, bool remove, int fromIndex, int toIndex) {
    if (remove) {
        _rM.rmvEdge(fromIndex, toIndex);
        if (!directed)
            _rM.rmvEdge(toIndex, fromIndex);
    } else {
        _rM.addEdge(fromIndex, toIndex);
        if (!directed)
            _rM.addEdge(toIndex, fromIndex);
    }
    _edgeFTInited = false;
    _fTEdgeInited = false;
}

int
convert2MiniZinc::getEdgeIDByIndex(int fromIndex, int toIndex) {
    int nR = _regionIDsByName.size();
    if (!_fTEdgeInited) {
        int id = 0;
        // update the eid to from-to pairs
        _fromToEdge.resize(nR, std::vector<int>(nR, -1));
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                _fromToEdge[j][outEdges[k]] = id;
                id++;
            }
        }
        _fTEdgeInited = true;
    }
    if ((fromIndex < nR) && (fromIndex >= 0) && (toIndex < nR) && (toIndex >= 0)) {
        return _fromToEdge[fromIndex][toIndex];
    } else {
        std::cout << "Query for wrong region ID.  Return meaningless data" << std::endl;
        return -1;
    }
}

int
convert2MiniZinc::getEdgeIDByName(std::string nameFrom, std::string nameTo) {
    return getEdgeIDByIndex(_regionIDsByName[nameFrom], _regionIDsByName[nameTo]);
}

std::pair<int, int>
convert2MiniZinc::regionFromToByEdgeIndex(int eid) {
    if (!_edgeFTInited) {
        int nR = _regionIDsByName.size();
        // update the eid to from-to pairs
        _edgeFromTo.clear();
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::pair<int, int> eFT = std::make_pair(j, outEdges[k]);
                _edgeFromTo.push_back(eFT);
            }
        }
        _edgeFTInited = true;
    }
    if (eid < _edgeFromTo.size()) {
        return _edgeFromTo[eid];
    } else {
        std::cout << "Queue for wrong edge ID.  Return meaningless data" << std::endl;
        return std::make_pair(-1, -1);
    }
}

void
convert2MiniZinc::printMiniZinc(printDest dest, printFor reason, std::string fileName, unsigned int layerNumber, std::vector<std::vector<int>> assignment, bool finalBC) {
    std::string miniZinc;
    switch (reason) {
    case (patching):
        if (layerNumber == 0)
            return;
        miniZinc = miniZincString(layerNumber, finalBC);
        break;

    case (reassignment):
        miniZinc = miniZincReassign(assignment, finalBC);
        break;

    case (planassignment):
        if (assignment.size() == 1)
            // will be a problem when the assignment has no transitions
            return;
        miniZinc = miniZincAssign(assignment, finalBC);
        break;

    default:
        std::cout << "The printFor is not implemented yet" << std::endl;
        break;
    }
    std::ofstream out(fileName + ".mzn");
    // for redirect the std out of the child process in the toFork case
    int link[2];
    // get the input/output file names done
    int fl = fileName.size();
    char mznFileName[fl + 5];
    char outFileName[fl + 5];
    char mznExecName[9];
    sprintf(mznFileName, "%s.mzn", fileName.c_str());
    sprintf(outFileName, "%s.txt", fileName.c_str());
    sprintf(mznExecName, "minizinc");
    // check where the mzn is printed to
    switch (dest) {
    case toTerminal:
        std::cout << "Printing " << fileName << " to Terminal: " << std::endl;
        std::cout << miniZinc << std::endl;
        break;

    case toFile:
        out << miniZinc;
        out.close();
        break;

    case toFork: {
        // print to file
        out << miniZinc;
        out.close();
        // call Minizinc
        link[1] = open(outFileName, O_WRONLY | O_CREAT, 0777);   // write side of the pipe, if it is used
        pid_t pid = fork();
        if (pid == 0) {
            // then it is the child
            char *argv[3] = {mznExecName, mznFileName, NULL};
            // re-direct the std out
            // std::cout << mznFileName << ", " << outFileName << std::endl;
            int resd = dup2(link[1], STDOUT_FILENO);
            execv("/home/brad/Downloads/MiniZincIDE-2.7.6-bundle-linux-x86_64/bin/minizinc", argv);
            // execvp("minizinc", argv);
        }
        // wait for the child to complete
        int res = waitPidHandy(pid);
        close(link[1]);
        break;
    }

    case toPatcher: {
        // print to file
        out << miniZinc;
        out.close();
        // call Minizinc
        // link[1] = memfd_create("read2Patch", 0);
        if (pipe(link) == -1) {
            std::cout << "Pipe creation failed" << std::endl;
            break;
        }
        std::cout << "to patcher" << std::endl;
        // FILE *fp = fdopen(link[1], O_RDONLY);
        pid_t pid1 = fork();
        if (pid1 == 0) {
            // then it is the child
            char *argv[3] = {mznExecName, mznFileName, NULL};
            // re-direct the std out
            // std::cout << mznFileName << ", " << outFileName << std::endl;
            int resd = dup2(link[1], STDOUT_FILENO);
            execv("/home/brad/Downloads/MiniZincIDE-2.7.6-bundle-linux-x86_64/bin/minizinc", argv);
            // execvp("minizinc", argv);
        }
        // get the string from pipe
        int bufSize = 100, ret = 0;
        char buff[bufSize + 1];
        buff[bufSize] = '\0';
        std::string mnzEop = "----------";
        // int ret = read(link[0], buff, bufSize);
        std::stringstream bufStream;
        do {
            ret = read(link[0], buff, bufSize);
            buff[ret] = '\0';
            bufStream << buff;
        } while (ret == bufSize);
        // std::cout << bufStream.str();
        // wait for the child to complete
        int res = waitPidHandy(pid1);
        close(link[1]);
        close(link[0]);
        // transform the string into _patch
        _patch.clear();
        std::string line;
        while (getline(bufStream, line)) {
            if (line.compare("=====UNSATISFIABLE=====") == 0) {
                break;
            }
            if (line.compare("----------") == 0) {
                break;
            }
            std::stringstream ss(line);
            std::string step, eN, brac;
            int sbt, sat, eA;
            ss >> step >> sbt >> sat >> brac;
            std::vector<int> state;
            while (!ss.eof()) {
                ss >> eN >> eA;
                state.push_back(eA);
            }
            _patch.push_back(state);
        }
        // print to debug
        std::cout << "Patch size: " << _patch.size() << std::endl;
        for (int i = 0; i < _patch.size(); i++) {
            for (int j = 0; j < _patch[i].size(); j++) {
                std::cout << _patch[i][j] << " ";
            }
            std::cout << std::endl;
        }
        break;
    }

    default:
        std::cout << "Print Destination " << dest << " not implemented yet" << std::endl;
        break;
    }
}

// utilities in printing miniZinc
// comment
std::string
mZCom(std::string comment) {
    std::ostringstream ret;
    ret << "% " << comment << std::endl;
    return ret.str();
}

// variable with range
std::string
mZVar(std::string name, int low, int high) {
    std::ostringstream ret;
    ret << "var " << low << ".." << high << ": " << name << "; " << std::endl;
    return ret.str();
}

// constraint
std::string
mZCtr(std::string constraint) {
    std::ostringstream ret;
    ret << "constraint " << constraint << "; " << std::endl;
    return ret.str();
}

// literals: true -> (X > 0); false -> (X == 0)
std::string
mznLiteral(bool n, std::string X) {
    if (n) {
        return "(" + X + " > 0)";
    } else {
        return "(" + X + " == 0)";
    }
}

std::string
mznEdgeLiteral(bool n, std::string X1, std::string X2) {
    if (n) {
        return "(" + X1 + " > 0) \\/ (" + X2 + " > 0)";
    } else {
        return "((" + X1 + " == 0) /\\ (" + X2 + " == 0))";
    }
}

std::string
convert2MiniZinc::miniZincAssign(std::vector<std::vector<int>> assignment, bool finalBC) {
    // here the assignment is the region swarm state
    std::stringstream os;
    int nR = _regionNames.size();
    int nS = assignment.size();   // number of transitions
    int sum = 0;
    // get the total robot number in assignment
    for (int i = 0; i < nR; i++) {
        sum += _iniNs[i];
    }
    // print region variables
    os << mZCom("Declare Regional Variables");
    for (int i = 0; i < nS; i++) {
        for (int j = 0; j < nR; j++) {
            std::string rP = "r" + _regionNames[j] + std::to_string(i + 1);
            os << mZVar(rP, 0, _maxNs[j]);
        }
    }
    // print edge variables
    os << mZCom("Declare Edge Variables");
    for (int i = 0; i < nS - 1; i++) {
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                os << mZVar(rE, 0, std::min(_maxNs[j], _maxNs[outEdges[k]]));
            }
        }
    }
    // Use EXC so that the problem is formulated as a optimization problem
    // os << mZVar("EXC", 0, sum);
    // print Initial/Final Condition
    os << mZCom("Initial-Final Conditions");
    for (int i = 0; i < nR; i++) {
        // initial region state = initial values
        std::string rPi = "r" + _regionNames[i] + std::to_string(1);
        os << mZCtr(rPi + " == " + std::to_string(_iniNs[i]));
        if (finalBC) {
            std::string rPf = "r" + _regionNames[i] + std::to_string(nS);
            os << mZCtr(rPf + " == " + std::to_string(_fnlNs[i]));
        }
    }
    // print transitions based on assignment and conservation rules
    os << mZCom("Conservation Law");
    for (int t = 0; t < nS; t++) {
        for (int i = 0; i < nR; i++) {
            // sum of out edges = sum of in edges for all regions
            if (t != (nS - 1)) {
                std::string si;
                std::vector<int> outNgh = getRegionOutNeighborByIndex(i);
                for (int j = 0; j < outNgh.size(); j++) {
                    std::string rE = _regionNames[i] + std::to_string(t + 1) + _regionNames[outNgh[j]] + std::to_string(t + 2);
                    si += rE + " + ";
                }
                std::string rP = "r" + _regionNames[i] + std::to_string(t + 1);
                os << mZCtr(si + " 0 == " + rP);
            }
            if (t != 0) {
                std::string si;
                std::vector<int> inNgh = getRegionInNeighborByIndex(i);
                for (int j = 0; j < inNgh.size(); j++) {
                    std::string rE = _regionNames[inNgh[j]] + std::to_string(t) + _regionNames[i] + std::to_string(t + 1);
                    si += rE + " + ";
                }
                std::string rP = "r" + _regionNames[i] + std::to_string(t + 1);
                os << mZCtr(si + " 0 == " + rP);
            }
        }
    }
    // print the region that is allowed or not allowed to have robot in each state
    os << mZCom("Symbolic Plan");
    for (int s = 0; s < nS; s++) {
        std::string si;
        for (int i = 0; i < nR; i++) {
            std::string rP = "r" + _regionNames[i] + std::to_string(s + 1);
            if (assignment[s][i] == 1) {
                si += mznLiteral(1, rP);
            } else if (assignment[s][i] == 0) {
                si += mznLiteral(0, rP);
            } else {
                std::cout << "Unknown assignment status: " << assignment[s][i] << std::endl;
            }
            si += " /\\ ";
        }
        si = si.substr(0, si.size() - 3);   // take away the last "\/ "
        os << mZCtr(si);
    }
    // Goal Function: We want to minimize the assignment exceed the regional capacity
    // os << mZCom("Goal Function");
    // std::string si;
    // std::vector<int> eMax = getEdgeMax();
    // for (int i = 0; i < nT; i++) {
    //     for (int j = 0; j < nE; j++) {
    //         int eas = assignment[i][j];
    //         if (eas > 0) {
    //             std::pair<int, int> ft = regionFromToByEdgeIndex(j);
    //             std::string rE = _regionNames[ft.first] + std::to_string(i + 1) + _regionNames[ft.second] + std::to_string(i + 2);
    //             si += "max( 0, " + rE + " - " + std::to_string(eMax[j]) + " ) + ";
    //         }
    //     }
    // }
    // os << mZCtr("EXC == " + si + " 0");
    // print the "Solve" and "Print"
    // os << "solve minimize EXC;" << std::endl;
    os << "solve satisfy;" << std::endl;
    for (int i = 0; i < nS - 1; i++) {
        os << "output [\"[step " << i + 1 << " " << i + 2 << " ] ";
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                os << _regionNames[j] << _regionNames[outEdges[k]] << ": \\(" << rE << ") ";
            }
        }
        os << "\\n\"]; " << std::endl;
    }
    // return
    return os.str();
}

std::string
convert2MiniZinc::miniZincReassign(std::vector<std::vector<int>> assignment, bool finalBC) {
    // here the assignment is the edge robot number assignment
    std::stringstream os;
    int nR = _regionNames.size();
    int nT = assignment.size();      // number of transitions
    int nE = assignment[0].size();   // number of edges in the road map
    int sum = 0;
    // get the total robot number in assignment
    for (int i = 0; i < nE; i++) {
        sum += assignment[0][i];
    }
    // print transition variables
    os << mZCom("Declare Transitions Variable");
    for (int i = 0; i < nT; i++) {
        for (int j = 0; j < nE; j++) {
            int eas = assignment[i][j];
            if (eas > 0) {
                std::pair<int, int> ft = regionFromToByEdgeIndex(j);
                std::string rE = _regionNames[ft.first] + std::to_string(i + 1) + _regionNames[ft.second] + std::to_string(i + 2);
                os << mZVar(rE, 1, sum);
            }
        }
    }
    os << mZVar("EXC", 0, sum);
    // print Initial/Final Condition
    os << mZCom("Initial-Final Conditions");
    for (int i = 0; i < nR; i++) {
        // initial region state = sum of out edges
        std::string si;
        std::vector<int> outNgh = getRegionOutNeighborByIndex(i);
        for (int j = 0; j < outNgh.size(); j++) {
            int eid = getEdgeIDByIndex(i, outNgh[j]);
            if (assignment[0][eid] > 0) {
                std::string rE = _regionNames[i] + std::to_string(1) + _regionNames[outNgh[j]] + std::to_string(2);
                si += rE + " + ";
            }
        }
        os << mZCtr(si + " 0 == " + std::to_string(_iniNs[i]));
        if (finalBC) {
            // final region state = sum of in edges
            // no need in the non - periodic case
            si = "";
            std::vector<int> inNgh = getRegionInNeighborByIndex(i);
            for (int j = 0; j < inNgh.size(); j++) {
                int eid = getEdgeIDByIndex(inNgh[j], i);
                if (assignment[nT - 1][eid] > 0) {
                    std::string rE = _regionNames[inNgh[j]] + std::to_string(nT) + _regionNames[i] + std::to_string(nT + 1);
                    si += rE + " + ";
                }
            }
            os << mZCtr(si + " 0 == " + std::to_string(_fnlNs[i]));
        }
    }
    // print transitions based on assignment and conservation rules
    os << mZCom("Conservation Law");
    for (int t = 1; t < nT; t++) {
        for (int i = 0; i < nR; i++) {
            // sum of out edges = sum of in edges for all regions
            std::string si;
            std::vector<int> outNgh = getRegionOutNeighborByIndex(i);
            for (int j = 0; j < outNgh.size(); j++) {
                int eid = getEdgeIDByIndex(i, outNgh[j]);
                if (assignment[t][eid] > 0) {
                    std::string rE = _regionNames[i] + std::to_string(t + 1) + _regionNames[outNgh[j]] + std::to_string(t + 2);
                    si += rE + " + ";
                }
            }
            si += " 0 == ";
            std::vector<int> inNgh = getRegionInNeighborByIndex(i);
            for (int j = 0; j < inNgh.size(); j++) {
                int eid = getEdgeIDByIndex(inNgh[j], i);
                if (assignment[t - 1][eid] > 0) {
                    std::string rE = _regionNames[inNgh[j]] + std::to_string(t) + _regionNames[i] + std::to_string(t + 1);
                    si += rE + " + ";
                }
            }
            os << mZCtr(si + " 0");
        }
    }
    // Goal Function: We want to minimize the assignment exceed the regional capacity
    os << mZCom("Goal Function");
    std::string si;
    std::vector<int> eMax = getEdgeMax();
    for (int i = 0; i < nT; i++) {
        for (int j = 0; j < nE; j++) {
            int eas = assignment[i][j];
            if (eas > 0) {
                std::pair<int, int> ft = regionFromToByEdgeIndex(j);
                std::string rE = _regionNames[ft.first] + std::to_string(i + 1) + _regionNames[ft.second] + std::to_string(i + 2);
                si += "max( 0, " + rE + " - " + std::to_string(eMax[j]) + " ) + ";
            }
        }
    }
    os << mZCtr("EXC == " + si + " 0");
    // print the "Solve" and "Print"
    os << "solve minimize EXC;" << std::endl;
    for (int i = 0; i < nT; i++) {
        os << "output [\"[step " << i + 1 << " " << i + 2 << " ] ";
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                int eid = getEdgeIDByIndex(j, outEdges[k]);
                int eas = assignment[i][eid];
                if (eas > 0) {
                    os << _regionNames[j] << _regionNames[outEdges[k]] << ": \\(" << rE << ") ";
                } else {
                    os << _regionNames[j] << _regionNames[outEdges[k]] << ": 0";
                }
            }
        }
        os << "\\n\"]; " << std::endl;
    }
    // return
    return os.str();
}

std::string
convert2MiniZinc::miniZincString(unsigned int layernumber, bool finalBC) {
    std::ostringstream os;
    // print region variables
    os << mZCom("Declare Regional Variables");
    int nR = _regionIDsByName.size();
    for (int i = 0; i <= layernumber; i++) {
        for (int j = 0; j < nR; j++) {
            std::string rP = "r" + _regionNames[j] + std::to_string(i + 1);
            os << mZVar(rP, 0, _maxNs[j]);
        }
    }
    // print edge variables
    os << mZCom("Declare Edge Variables");
    for (int i = 0; i < layernumber; i++) {
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                os << mZVar(rE, 0, std::min(_maxNs[j], _maxNs[outEdges[k]]));
            }
        }
    }
    // print clauses
    os << mZCom("Strategy Constraints");
    int nC = _clauses.size();
    for (int i = 0; i < layernumber; i++) {
        for (int j = 0; j < nC; j++) {
            int jB = _clauses[j].beforeRegionIndex.size();
            std::string cj;
            for (int k = 0; k < jB; k++) {
                int iR = _clauses[j].beforeRegionIndex[k];
                bool cm = _clauses[j].beforeComplemented[k];
                std::string rP = "r" + _regionNames[iR] + std::to_string(i + 1);
                std::string li = mznLiteral(!cm, rP);
                cj += li + " \\/ ";
            }
            int jA = _clauses[j].afterRegionIndex.size();
            for (int k = 0; k < jA; k++) {
                int iR = _clauses[j].afterRegionIndex[k];
                bool cm = _clauses[j].afterComplemented[k];
                std::string rP = "r" + _regionNames[iR] + std::to_string(i + 2);
                std::string li = mznLiteral(!cm, rP);
                cj += li + " \\/ ";
            }
            cj = cj.substr(0, cj.size() - 3);   // take away the last "\/ "
            os << mZCtr(cj);
        }
    }
    // print intermediate state literals
    /*
    os << mZCom("Intermediate State Constraints");
    nC = _edgeLiterals.size();   // number of clauses
    for (int i = 0; i < layernumber; i++) {
        for (int j = 0; j < nC; j++) {
            int jC = _edgeLiterals[j].size();   // number of literal in the jth clauses
            std::string cj;
            for (int k = 0; k < jC; k++) {
                std::pair<int, bool> ljk = _edgeLiterals[j][k];   // the literal
                std::pair<int, int> eid = regionFromToByEdgeIndex(ljk.first);
                std::string eP1 = _regionNames[eid.first] + std::to_string(i + 1) + _regionNames[eid.second] + std::to_string(i + 2);
                std::string eP2 = _regionNames[eid.second] + std::to_string(i + 1) + _regionNames[eid.first] + std::to_string(i + 2);
                std::string li = mznEdgeLiteral(ljk.second, eP1, eP2);
                cj += li + " \\/ ";
            }
            cj = cj.substr(0, cj.size() - 3);   // take away the last "\/ "
            os << mZCtr(cj);
        }
    }
    */
    // print conservation law: inlet sum equals after region value, and outlet sum equals before region value
    os << mZCom("Conservation Constraint: Inlets");
    for (int i = 0; i < layernumber; i++) {
        for (int j = 0; j < nR; j++) {
            std::string sj;
            auto inEdges = _rM.getInRegion(j);
            for (int k = 0; k < inEdges.size(); k++) {
                std::string rE = _regionNames[inEdges[k]] + std::to_string(i + 1) + _regionNames[j] + std::to_string(i + 2);
                sj += (rE + " + ");
            }
            // check if the inlet regions exist
            if (sj.size() > 2) {
                sj = sj.substr(0, sj.size() - 2);   // take away the last "+ "
                std::string rP = "r" + _regionNames[j] + std::to_string(i + 2);
                os << mZCtr(sj + "== " + rP);
            }
        }
    }
    os << mZCom("Conservation Constraint: Outlets");
    for (int i = 0; i < layernumber; i++) {
        for (int j = 0; j < nR; j++) {
            std::string sj;
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                sj += (rE + " + ");
            }
            // check if the outlet regions exist
            if (sj.size() > 2) {
                sj = sj.substr(0, sj.size() - 2);   // take away the last "+ "
                std::string rP = "r" + _regionNames[j] + std::to_string(i + 1);
                os << mZCtr(sj + "== " + rP);
            }
        }
    }
    // print initial and final condition
    os << mZCom("Initial Conditions");
    for (int j = 0; j < nR; j++) {
        std::string rP = "r" + _regionNames[j] + "1";
        os << mZCtr(rP + " == " + std::to_string(_iniNs[j]));
    }
    if (finalBC) {
        os << mZCom("Final Conditions");
        for (int j = 0; j < nR; j++) {
            std::string rP = "r" + _regionNames[j] + std::to_string(layernumber + 1);
            os << mZCtr(rP + " == " + std::to_string(_fnlNs[j]));
        }
    }
    // print "solve" and output
    os << "solve satisfy;" << std::endl;
    for (int i = 0; i < layernumber; i++) {
        os << "output [\"[step " << i + 1 << " " << i + 2 << " ] ";
        for (int j = 0; j < nR; j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                std::string rE = _regionNames[j] + std::to_string(i + 1) + _regionNames[outEdges[k]] + std::to_string(i + 2);
                os << _regionNames[j] << _regionNames[outEdges[k]] << ": \\(" << rE << ") ";
            }
        }
        os << "\\n\"]; " << std::endl;
    }
    // return
    return os.str();
}

std::string
convert2MiniZinc::regionState2NodeString(std::vector<int> rS, int stateNum) {
    std::string ret;
    if (rS.size() == 0) {
        std::cout << "State vector size is 0.  Problems in regionState2NodeString.  " << std::endl;
        return ret;
    }
    for (int i = 0; i < rS.size(); i++) {
        ret += _regionNames[i] + std::to_string(stateNum) + " [label=\"" + _regionNames[i] + ": " + std::to_string(rS[i]) + "\"];\r\n";
    }
    // specify options so that the node is at one rank and ordered via adding invisible edges
    ret += "rank = same {";
    bool first = true;
    for (int i = 0; i < rS.size(); i++) {
        if (first) {
            first = false;
            ret += _regionNames[i] + std::to_string(stateNum);
        } else {
            ret += " -> " + _regionNames[i] + std::to_string(stateNum);
        }
    }
    ret += " [style=invis] }";
    return ret;
}

void
convert2MiniZinc::printPatch2Dot(std::string fileName, std::vector<std::vector<int>> patch) {
    if (patch.size() == 0) {
        std::cout << "Patch size is 0.  Cancel printing to file.  " << std::endl;
        return;
    }
    std::ofstream out(fileName + ".dot");
    // print header and the first state
    out << "digraph {" << std::endl;
    std::vector<int> rS = edge2RegionState(true, patch[0]);
    out << regionState2NodeString(rS, 0);
    // for each transition
    for (int i = 0; i < patch.size(); i++) {
        // print the node in this state
        // find the number of assignment inside each region after this transition
        rS = edge2RegionState(false, patch[i]);
        out << regionState2NodeString(rS, i + 1);
        // print the edge connect to the last state
        int edgeID = 0;
        for (int j = 0; j < _regionNames.size(); j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                int fromJ2outK = patch[i][edgeID];   // before: sum of out edge
                if (fromJ2outK != 0) {
                    out << _regionNames[j] + std::to_string(i) + " -> " + _regionNames[outEdges[k]] + std::to_string(i + 1) << std::endl;
                }
                edgeID++;
            }
        }
    }
    // print invisible verticle skeleton
    for (int i = 0; i <= patch.size(); i++) {
        out << _regionNames[0] + std::to_string(i);
        if (i != patch.size()) {
            out << " -> ";
        }
    }
    out << " [ style=invis; weight=1000 ]" << std::endl;
    for (int i = 0; i <= patch.size(); i++) {
        out << _regionNames[_regionNames.size() - 1] + std::to_string(i);
        if (i != patch.size()) {
            out << " -> ";
        }
    }
    out << " [ style=invis; weight=1000 ]" << std::endl;
    // file ender
    out << "}" << std::endl;
}

void
convert2MiniZinc::printPatch2Dot(std::string fileName) {
    return printPatch2Dot(fileName, _patch);
}

void
convert2MiniZinc::printEdgeLiterals() {
    for (int i = 0; i < _edgeLiterals.size(); i++) {
        std::vector<std::pair<int, bool>> li = _edgeLiterals[i];
        std::cout << "[Clause " << i << "]: ";
        for (int j = 0; j < li.size(); j++) {
            std::pair<int, bool> pij = li[j];
            int eid;
            std::string From, To;
            eid = pij.first;
            std::pair<int, int> idFT = regionFromToByEdgeIndex(eid);
            From = _regionNames[idFT.first];
            To = _regionNames[idFT.second];
            if (pij.second) {
                std::cout << From << To << " \\/ ";
            } else {
                std::cout << " !" << From << To << " \\/ ";
            }
        }
        std::cout << std::endl;
    }
}

void
convert2MiniZinc::testPatchWithEdgeLiterals() {
    int nP = _patch.size();
    int nC = _edgeLiterals.size();
    for (int i = 0; i < nP; i++) {
        std::vector<int> pi = _patch[i];
        for (int j = 0; j < nC; j++) {
            std::vector<std::pair<int, bool>> cj = _edgeLiterals[j];
            bool same = false;
            for (int cji = 0; cji < cj.size(); cji++) {
                std::pair<int, bool> lit = cj[cji];
                std::pair<int, int> idFT = regionFromToByEdgeIndex(lit.first);
                int conjEid = getEdgeIDByIndex(idFT.second, idFT.first);
                if (lit.second) {
                    bool pred = false;
                    for (int pei = 0; pei < pi.size(); pei++) {
                        bool ei = (pi[pei] > 0);
                        if ((lit.first == pei) || (conjEid == pei)) {
                            pred |= (ei == lit.second);
                        }
                    }
                    same |= pred;
                } else {
                    bool pred = true;
                    for (int pei = 0; pei < pi.size(); pei++) {
                        bool ei = (pi[pei] > 0);
                        if ((lit.first == pei) || (conjEid == pei)) {
                            pred &= (ei == lit.second);
                        }
                    }
                    same |= pred;
                }
            }
            if (!same) {
                std::cout << "======================================" << std::endl;
                std::cout << "Conflict between transition and clause" << std::endl;
                std::cout << "[Transition" << i << "]: ";
                for (int pei = 0; pei < pi.size(); pei++) {
                    std::pair<int, int> idFT = regionFromToByEdgeIndex(pei);
                    std::string From, To;
                    From = _regionNames[idFT.first];
                    To = _regionNames[idFT.second];
                    std::cout << From << " -> " << To << ": " << pi[pei] << ", ";
                }
                std::cout << std::endl;
                std::cout << "[Clause " << j << "]: ";
                for (int cji = 0; cji < cj.size(); cji++) {
                    std::pair<int, bool> pij = cj[cji];
                    int eid;
                    std::string From, To;
                    eid = pij.first;
                    std::pair<int, int> idFT = regionFromToByEdgeIndex(eid);
                    From = _regionNames[idFT.first];
                    To = _regionNames[idFT.second];
                    if (pij.second) {
                        std::cout << From << To << " \\/ ";
                    } else {
                        std::cout << " !" << From << To << " \\/ ";
                    }
                }
                std::cout << std::endl;
            }
        }
    }
}