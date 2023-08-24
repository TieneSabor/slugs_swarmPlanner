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
}

void
convert2MiniZinc::printMiniZinc(printDest dest, std::string fileName, unsigned int layerNumber) {
    if (layerNumber == 0)
        return;
    std::string miniZinc = miniZincString(layerNumber);
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
        std::cout << "to patcher" << std::endl;
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
        std::cout << "to patcher" << std::endl;
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
convert2MiniZinc::miniZincString(unsigned int layernumber) {
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
    os << mZCom("Final Conditions");
    for (int j = 0; j < nR; j++) {
        std::string rP = "r" + _regionNames[j] + std::to_string(layernumber + 1);
        os << mZCtr(rP + " == " + std::to_string(_fnlNs[j]));
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
convert2MiniZinc::printPatch2Dot(std::string fileName) {
    if (_patch.size() == 0) {
        std::cout << "Patch size is 0.  Cancel printing to file.  " << std::endl;
        return;
    }
    std::ofstream out(fileName + ".dot");
    // print header and the first state
    out << "digraph {" << std::endl;
    std::vector<int> rS = edge2RegionState(true, _patch[0]);
    out << regionState2NodeString(rS, 0);
    // for each transition
    for (int i = 0; i < _patch.size(); i++) {
        // print the node in this state
        // find the number of assignment inside each region after this transition
        rS = edge2RegionState(false, _patch[i]);
        out << regionState2NodeString(rS, i + 1);
        // print the edge connect to the last state
        int edgeID = 0;
        for (int j = 0; j < _regionNames.size(); j++) {
            auto outEdges = _rM.getOutRegion(j);
            for (int k = 0; k < outEdges.size(); k++) {
                int fromJ2outK = _patch[i][edgeID];   // before: sum of out edge
                if (fromJ2outK != 0) {
                    out << _regionNames[j] + std::to_string(i) + " -> " + _regionNames[outEdges[k]] + std::to_string(i + 1) << std::endl;
                }
                edgeID++;
            }
        }
    }
    // print invisible verticle skeleton
    for (int i = 0; i <= _patch.size(); i++) {
        out << _regionNames[0] + std::to_string(i);
        if (i != _patch.size()) {
            out << " -> ";
        }
    }
    out << " [ style=invis; weight=1000 ]" << std::endl;
    for (int i = 0; i <= _patch.size(); i++) {
        out << _regionNames[_regionNames.size() - 1] + std::to_string(i);
        if (i != _patch.size()) {
            out << " -> ";
        }
    }
    out << " [ style=invis; weight=1000 ]" << std::endl;
    // file ender
    out << "}" << std::endl;
}