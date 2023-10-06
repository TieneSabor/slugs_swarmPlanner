#include <extensionSwarmTest.hpp>

using namespace swarmTest;

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::init(std::list<std::string> &filenames) {
    // build the regional map
    BF RMSG = parseRegionMap(filenames);
    // parse all GR1 specs except for the edge safety guarantee
    T::init(filenames);
    if (filenames.size() == 0) {
        std::cerr << "Error: Need a file name for extracting a symbolic "
                     "strategy.\n";
        throw "Please adapt the parameters.";
    } else {
        outputFilename = filenames.front();
        filenames.pop_front();
    }
    // add region map safety guarantee to safetySys
    safetySysNoRM = safetySys;
    safetySys &= RMSG;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::getVarFromCVZID(bool Pre, int cvzID) {
    if (Pre) {
        return variables[2 * cvzID];
    } else {
        return variables[2 * cvzID + 1];
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::regAsn2BF(std::vector<int> rA) {
    BF ret = mgr.constantTrue();
    // slug pre-var  id = 2*RM id in cvz
    for (int i = 0; i < pcvz->getRegionNum(); i++) {
        if (rA[i] > 0) {
            ret &= variables[2 * i];
        } else {
            ret &= !variables[2 * i];
        }
    }
    return ret;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::getRMSG() {
    // slug pre-var  id = 2*RM id in cvz
    // slug post-var id = 2*RM id in cvz + 1
    // add edges as safety guarantee
    BF RMSG = mgr.constantTrue();
    for (int i = 0; i < pcvz->getRegionNum(); i++) {
        std::vector<int> inIds = pcvz->getRegionInNeighborByIndex(i);
        BF inSG = !variables[2 * i + 1];
        for (int j = 0; j < inIds.size(); j++) {
            inSG |= variables[2 * inIds[j]];
        }
        RMSG &= inSG;
        std::vector<int> outIds = pcvz->getRegionOutNeighborByIndex(i);
        BF outSG = !variables[2 * i];
        for (int j = 0; j < outIds.size(); j++) {
            outSG |= variables[2 * outIds[j] + 1];
        }
        RMSG &= outSG;
    }
    return RMSG;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::parseRegionMap(std::list<std::string> &filenames) {
    std::string inFileName = filenames.front();
    // filenames.pop_front();
    // Open input file or produce error message if that does not work
    std::ifstream inFile(inFileName.c_str());
    int readMode = -1;
    std::string currentLine;
    std::vector<std::string> regs;
    std::vector<int> caps;
    std::vector<int> inis;
    std::vector<int> fnls;
    while (std::getline(inFile, currentLine)) {
        if (0 == (currentLine.compare("# specify regions"))) {
            readMode = 1;
        } else if (0 == (currentLine.compare("# specify edges"))) {
            // initialize cvz class
            pcvz = new convert2MiniZinc(regs);
            for (int i = 0; i < regs.size(); i++) {
                pcvz->setRegionConstByName(regs[i], caps[i], inis[i], fnls[i]);
            }
            readMode = 2;
        } else if (0 == (currentLine.compare("# specify decreased region max number"))) {
            readMode = 3;
        } else if (0 == (currentLine.compare("# specify removed edge"))) {
            readMode = 4;
        } else if (0 == (currentLine.compare("# reassign subswarm number"))) {
            readMode = 5;
        } else if (0 == (currentLine.compare("# regional map end"))) {
            readMode = 0;
            break;
        } else {
            std::istringstream is(currentLine);
            std::string shrp, reg, cap, ini, fnl, from, to, num, reasgnSN, reasgnGN;
            switch (readMode) {
            case (1):
                // get the region name, capacity, initial number and final number
                is >> shrp >> reg >> cap >> ini >> fnl;
                regs.push_back(reg);
                caps.push_back(stoi(cap));
                inis.push_back(stoi(ini));
                fnls.push_back(stoi(fnl));
                break;
            case (2):
                is >> shrp >> from >> to;
                if (0 == from.compare(to)) {
                    // self loop
                    pcvz->addEdgeByName(true, false, from, to);
                } else {
                    // assume all inter-region edges are bi-directional
                    pcvz->addEdgeByName(false, false, from, to);
                }
                break;
            case (3): {
                is >> shrp >> reg >> cap;
                modReg = reg;
                modMax = stoi(cap);
                break;
            }
            case (4):
                is >> shrp >> from >> to;
                rmvEdgeFrom = from;
                rmvEdgeTo = to;
                break;
            case (5):
                is >> shrp >> reasgnSN >> reasgnGN;
                reasgnGoalNum = stoi(reasgnGN);
                reasgnStateNum = stoi(reasgnSN);
                for (int i = 0; i < regs.size(); i++) {
                    is >> num;
                    reasgnRegionState.push_back(stoi(num));
                }
                break;

            default:
                std::cout << "un defined regional map input mode: " << readMode << std::endl;
                break;
            }
        }
    }
    // add variables to slugs
    for (int i = 0; i < regs.size(); i++) {
        addVariable(PreOutput, regs[i]);
        addVariable(PostOutput, regs[i] + "'");
    }
    // slug pre-var  id = 2*RM id in cvz
    // slug post-var id = 2*RM id in cvz + 1
    initState = regAsn2BF(inis);
    return getRMSG();
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::updateLayers(int goalID, int layerNum, BF data) {
    int sl = layers2Goals[goalID].size();
    if (layerNum >= sl) {
        layers2Goals[goalID].resize(layerNum + 1);
        layers2Goals[goalID][layerNum] = data;
        return;
    }
    layers2Goals[goalID][layerNum] |= data;
    // dbg
    // std::cout << "Hi " << std::endl;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::updateL2G() {
    layers2Goals.resize(livenessGuarantees.size());
    // for (int j = 0; j < livenessGuarantees.size(); j++) {
    // }
    int livAs = 0, livGr = 0, level = 0;
    // livAs++ everytime
    // if livAs >= num of live. asum., reset livAs and level++
    // if next goal_number is different, reset level and livGr++
    // if livGr >= num of live. gura., reset livGr and go to next "Nu_2"
    bool nxtNu2 = false;
    for (int d = 0; d < strategyDumpingData.size(); d++) {
        // check "state machine of updateL2G"
        std::cout << "d: " << d << ", livAs: " << livAs << ", level: " << level << ", livGr: " << livGr << std::endl;
        // update layer
        BF state = strategyDumpingData[d].second.ExistAbstract(varCubePostOutput);   // no idea on pre/post input actually... Should not matter for now
        // dbg
        // std::cout << "Hi " << std::endl;
        updateLayers(livGr, level, state);
        livAs++;
        // update "state machine"
        if (livAs >= livenessAssumptions.size()) {
            livAs = 0;
            level++;
        }
        int nxtGoal = strategyDumpingData[d].first;
        if ((d + 1) < strategyDumpingData.size()) {
            nxtGoal = strategyDumpingData[d + 1].first;
        }
        if (nxtGoal != strategyDumpingData[d].first) {
            level = 0;
            livGr++;
        }
        if (livGr >= livenessGuarantees.size()) {
            livGr = 0;
        }
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
int
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::checkLayer(BF tbc, int goalID) {
    if (goalID >= layers2Goals.size()) {
        std::cout << "Goal ID not exists: " << goalID << std::endl;
        return -1;
    }
    for (int i = 0; i < layers2Goals[goalID].size(); i++) {
        if (!(layers2Goals[goalID][i] & tbc).isFalse()) {
            return i;
        }
    }
    return -1;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::printStrategyBDD(BF combinedStrategy, std::string filename) {
    std::ostringstream fileExtraHeader;
    fileExtraHeader << "# This file is a BDD exported by the SLUGS\n#\n# "
                       "This BDD is a strategy.\n";
    if (!systemGoalEncoded)
        fileExtraHeader << "#\n# This header contains extra information "
                           "used by LTLMoP's BDDStrategy.\n";
    fileExtraHeader << "# Currently, the only metadata is 1) the total "
                       "number of system goals\n";
    fileExtraHeader << "# and 2) the mapping between variable numbers and "
                       "proposition names.\n#\n";
    fileExtraHeader << "# Some special variables are also added:\n";
    fileExtraHeader << "#       - `_jx_b*` are used as a binary vector (b0 "
                       "is LSB) to indicate\n";
    fileExtraHeader
        << "#         the index of the currently-pursued goal.\n";
    if (!systemGoalEncoded) {
        fileExtraHeader << "#       - `strat_type` is a binary variable "
                           "used to indicate whether we are\n";
        fileExtraHeader << "#          moving closer to the current goal "
                           "(0) or transitioning to the next goal (1)\n#\n";
    }
    fileExtraHeader << "# Num goals: " << livenessGuarantees.size() << "\n";
    fileExtraHeader << "# Variable names:\n";
    if (!systemGoalEncoded) {
        for (unsigned int i = 0; i < variables.size(); i++) {
            fileExtraHeader << "#\t" << i << ": " << variableNames[i]
                            << "\n";
        }
    } else {
        for (unsigned int i = 0; i < variables.size(); i++) {
            if (doesVariableInheritType(i, PreInput)) {
                fileExtraHeader << "#\t" << i << ": in_" << variableNames[i]
                                << "\n";
            } else if (doesVariableInheritType(i, PostInput)) {
                fileExtraHeader << "#\t" << i << ": in_" << variableNames[i]
                                << "\n";
            } else {
                fileExtraHeader << "#\t" << i << ": " << variableNames[i]
                                << "\n";
            }
        }
    }
    fileExtraHeader
        << "#\n# For information about the DDDMP format, please see:\n";
    fileExtraHeader << "#    "
                       "http://www.cs.uleth.ca/~rice/cudd_docs/dddmp/"
                       "dddmpAllFile.html#dddmpDump.c\n#\n";
    fileExtraHeader << "# For information about how this file is "
                       "generated, please see the SLUGS source.\n#\n";

    // take away extra variables: strat, _jx_b*, ...
    // assign whether going to the next goal: no
    // combinedStrategy &= !variables[goalTransitionSelectorVar];
    // calculate encoding
    // int goalNumb = 0;
    // for (unsigned j = 0; j < counterVarNumbersPre.size(); j++) {
    //     if (goalNumb & (1 << j)) {
    //         combinedStrategy &= variables[counterVarNumbersPre[j]];
    //     } else {
    //         combinedStrategy &= !variables[counterVarNumbersPre[j]];
    //     }
    // }
    // do the exist abstract
    // ...
    // above section should produce the same thing as positionalStrategiesForTheIndividualGoals[goalNumb]
    // int goalNumb = 0;
    // combinedStrategy = positionalStrategiesForTheIndividualGoals[goalNumb];

#ifndef NDEBUG
    std::string tempFilename(tmpnam(NULL));
    tempFilename = tempFilename + "_strategyBdd.dot";
    std::cerr << "Writing DOT file of the BDD to: " << tempFilename
              << std::endl;
    BF_newDumpDot(*this, combinedStrategy,
                  "SymbolicStrategyCounterVar PreInput PreOutput PostInput "
                  "PostOutput",
                  tempFilename.c_str());
#endif

    mgr.writeBDDToFile(filename.c_str(), fileExtraHeader.str(),
                       combinedStrategy, variables, variableNames);
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::safetySys2IntStateRec(BF F, std::vector<std::pair<int, int>> edgePred, std::vector<std::pair<int, bool>> literal, std::vector<std::vector<std::pair<int, bool>>> &literals) {
    // std::cout << "level: " << edgePred.size() << std::endl;
    // Boundary Condition: If edgePred is empty, return
    if (F.isFalse()) {
        literals.push_back(literal);
        return;
    }
    if (edgePred.size() == 0) {
        // BF EAF = F.ExistAbstract(varCubePre);
        // printStrategyBDD(EAF, "remainBDD" + std::to_string(literals.size()));
        // if (!EAF.isTrue()) {
        //     std::cout << "Bad" << std::endl;
        // } else {
        //     // literals.push_back(literal);
        //     std::cout << "?" << std::endl;
        // }
        return;
    }
    // Check if the edgePred.back() existed in the regional map
    std::pair<int, int> e = edgePred.back();
    std::vector<std::pair<int, int>> nEP = edgePred;
    nEP.pop_back();
    int eid = pcvz->getEdgeIDByIndex(e.first, e.second);
    // If so, do the Boole Expansion: F = ("Left") (F(e=0) + e)*("Right") (F(e=1) + (!e))
    if (eid != -1) {
        std::vector<std::pair<int, bool>> leftLit, rightLit;
        BF leftF = (F & (!(variables[e.first * 2] & variables[e.second * 2])));
        BF rightF = F;   //(F & ((variables[e.first * 2] & variables[e.second * 2])));
        leftLit = literal;
        rightLit = literal;
        leftLit.push_back(std::make_pair(eid, false));
        rightLit.push_back(std::make_pair(eid, true));
        safetySys2IntStateRec(leftF, nEP, leftLit, literals);
        safetySys2IntStateRec(rightF, nEP, rightLit, literals);
    }
    // If not, pass: F = F(e=0)*(!e)
    else {
        BF nextF = (F & (!(variables[e.first * 2] & variables[e.second * 2])));
        safetySys2IntStateRec(nextF, nEP, literal, literals);
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::testingBF() {
    int idP = 5, idQ = 6, idR = 7, idS = 8;
    BF bfP = variables[idP * 2];
    BF bfQ = variables[idQ * 2];
    BF bfR = variables[idR * 2];
    BF bfS = variables[idS * 2];
    BF tbt = safetySysNoRM;
    BF sbf = (tbt & (!(bfS & bfS)));
    sbf = (sbf & (!(bfR & bfS)));
    sbf = (sbf & (!(bfR & bfR)));
    sbf = (sbf & ((bfQ & bfS)));
    if (sbf.isFalse()) {
        std::cout << "indeed" << std::endl;
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::safetySys2IntState() {
    // create intermediate state predicate variables in BF manager
    // Not required because for now we only care about CNF
    // Create a list of edges to process recursively, including those not in region map
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < pcvz->getRegionNum(); i++) {
        for (int j = i; j < pcvz->getRegionNum(); j++) {
            edges.push_back(std::make_pair(i, j));
        }
    }
    // transform safetySysNoRM to CNF in int. state
    std::vector<std::pair<int, bool>> literal;
    std::vector<std::vector<std::pair<int, bool>>> literals;
    BF tbt = !safetySysNoRM;   // to-be-transformed
    // printStrategyBDD(safetySysNoRM, "safetyNoRMBDD");
    safetySys2IntStateRec(tbt, edges, literal, literals);
    // Import CNF to pcvz
    pcvz->setEdgeLiterals(literals);
    // pcvz->printEdgeLiterals();
    // testingBF();
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::safetyTransitionSetRec(BF F, std::vector<std::pair<int, int>> edgePred) {
    // std::cout << "level: " << edgePred.size() << std::endl;
    // Boundary Condition: If edgePred is empty, return
    // if (F.isFalse()) {
    //     // BF or False = BF
    //     return mgr.constantFalse();
    // }
    if (edgePred.size() == 0) {
        // BF EAF = F.ExistAbstract(varCubePre);
        // printStrategyBDD(EAF, "remainBDD" + std::to_string(literals.size()));
        // if (!EAF.isTrue()) {
        //     std::cout << "Bad" << std::endl;
        // } else {
        //     // literals.push_back(literal);
        //     std::cout << "?" << std::endl;
        // }
        // This will make the clause vanish: BF1 and (BF2 or True) = BF1 and True = BF1
        // return mgr.constantTrue();
        if (F.isFalse()) {
            // BF or False = BF
            return mgr.constantFalse();
        } else {
            // This will make the clause vanish: BF1 and (BF2 or True) = BF1 and True = BF1
            return mgr.constantTrue();
        }
    }
    // Check if the edgePred.back() existed in the regional map
    std::pair<int, int> e = edgePred.back();
    std::vector<std::pair<int, int>> nEP = edgePred;
    nEP.pop_back();
    int eid = pcvz->getEdgeIDByIndex(e.first, e.second);
    BF EF = (variables[e.first * 2] & variables[e.second * 2 + 1]) | (variables[e.first * 2 + 1] & variables[e.second * 2]);
    // If so, do the Boole Expansion: F = ("Left") (F(e=0) + e)*("Right") (F(e=1) + (!e))
    if (eid != -1) {
        std::vector<std::pair<int, bool>> leftLit, rightLit;
        BF leftF = F;   // (F & (!(variables[e.first * 2] & variables[e.second * 2])));
        BF rightF = (F & ((variables[e.first * 2] & variables[e.second * 2])));
        BF leftAns = safetyTransitionSetRec(leftF, nEP);
        BF rightAns = safetyTransitionSetRec(rightF, nEP);
        return (leftAns | (EF)) & (rightAns | (!EF));
    }
    // If not, pass: F = F(e=0)*(!e)
    else {
        BF nextF = F;   // (F & (!(variables[e.first * 2] & variables[e.second * 2])));
        BF ans = safetyTransitionSetRec(nextF, nEP);
        return ans;   // & (!EF);
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::safetyTransitionSet() {
    // create intermediate state predicate variables in BF manager
    // Not required because for now we only care about CNF
    // Create a list of edges to process recursively, including those not in region map
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i < pcvz->getRegionNum(); i++) {
        for (int j = i; j < pcvz->getRegionNum(); j++) {
            edges.push_back(std::make_pair(i, j));
        }
    }
    // transform safetySysNoRM to CNF in int. state
    BF tbt = safetySysNoRM;   // to-be-transformed
    // printStrategyBDD(safetySysNoRM, "safetyNoRMBDD");
    BF ret = safetyTransitionSetRec(tbt, edges);
    // verify intermediate state checker
    // for (int i = 0; i < edges.size(); i++) {
    //     for (int j = i + 1; j < edges.size(); j++) {
    //         auto ei = edges[i];
    //         auto ej = edges[j];
    //         BF EFi = (variables[ei.first * 2] & variables[ei.second * 2 + 1]) | (variables[ei.first * 2 + 1] & variables[ei.second * 2]);
    //         BF EFj = (variables[ej.first * 2] & variables[ej.second * 2 + 1]) | (variables[ej.first * 2 + 1] & variables[ej.second * 2]);
    //         BF check = EFi & EFj & ret;
    //         int eidi = pcvz->getEdgeIDByIndex(ei.first, ei.second);
    //         int eidj = pcvz->getEdgeIDByIndex(ej.first, ej.second);
    //         if ((eidj != -1) && (eidi != -1) && (check.isFalse())) {
    //             std::cout << "No edge " << variableNames[ei.first * 2] << " - " << variableNames[ei.second * 2] << " and "
    //                       << variableNames[ej.first * 2] << " - " << variableNames[ej.second * 2] << std::endl;
    //         }
    //     }
    // }
    return ret;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
BF
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::getIntermediateStateFeedback(std::vector<std::vector<int>> transitions) {
    BF fb = mgr.constantTrue();
    // iterate for each transitions
    int nT = transitions.size();
    for (int t = 0; t < nT; t++) {
        std::vector<int> lastState = pcvz->edge2RegionState(true, transitions[t]);
        std::vector<int> thisState = pcvz->edge2RegionState(false, transitions[t]);
        // check for intermediate state
        // phi_tau_1
        BF phi_tau_1 = mgr.constantTrue();
        for (int li = 0; li < lastState.size(); li++) {
            for (int ti = 0; ti < thisState.size(); ti++) {
                int liAs = lastState[li], tiAs = thisState[ti];
                if ((liAs > 0) && (tiAs > 0)) {
                    phi_tau_1 &= getVarFromCVZID(true, li) | getVarFromCVZID(true, ti);
                }
            }
        }
        // phi_tau_2
        BF phi_tau_2 = mgr.constantTrue();
        for (int si = 0; si < lastState.size(); si++) {
            if ((lastState[si] == 0) && (thisState[si] == 0)) {
                phi_tau_2 &= !getVarFromCVZID(true, si);
            }
        }
        // phi_q_i
        BF phi_q_i = regAsn2BF(lastState);
        BF phi_q_iP1 = regAsn2BF(thisState);
        // check that
        BF checker = phi_tau_1 & phi_tau_2 & (!phi_q_i) & (!phi_q_iP1);
        // if not passed, formulate feedbacks
        if (!(checker & (!safetySysNoRM)).isFalse()) {
            std::cout << "Transition " << t << " has dangerous intermediate state." << std::endl;
            pcvz->printTransition(transitions[t]);
            BF transBF = phi_q_i & (phi_q_iP1.SwapVariables(varVectorPre, varVectorPost));
            fb &= !transBF;
        }
    }
    return fb;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::BF2CVZClauses(BF bf, whichClause wc) {
    std::string rafilename = "bf_2_cnf";
    std::string header_ra_cnf = "";
    // TODO: consider making clause buffer a global thing
    int rclausesInCNFIndex[1000000], rclauseN = 0, rvarSize = 0;
    mgr.writeBDDToCNFFile(rafilename.c_str(), header_ra_cnf,
                          bf, variables, variableNames, rclausesInCNFIndex, rclauseN, rvarSize);
    for (int i = 0; i < rclauseN; i++) {
        pcvz->newClause(wc);
        // std::cout << "Clause #" << i << " ";
        for (int j = 0; j < rvarSize; j++) {
            // std::cout << clausesInCNFIndex[i * varSize + j] << " ";
            std::string lit = variableNames[j];
            switch (rclausesInCNFIndex[i * rvarSize + j]) {
            case 0:
                pcvz->addLiteral2LastClauseByName(wc, lit);
                break;
            case 1:
                pcvz->addLiteral2LastClauseByName(wc, "-" + lit);
                break;
            case 2:
                break;
            default:
                std::cout << "Unknown literal value: " << rclausesInCNFIndex[i * rvarSize + j] << std::endl;
            }
        }
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
std::vector<BF>
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::computeSymbolicStrategy() {
    // Prepare initial to-do list from the allowed initial states
    BF init = (oneStepRecovery) ? (winningPositions & initSys)
                                : (winningPositions & initSys & initEnv);

    // Prepare positional strategies for the individual goals
    std::vector<BF> positionalStrategiesForTheIndividualGoals(
        livenessGuarantees.size());
    for (unsigned int i = 0; i < livenessGuarantees.size(); i++) {
        BF casesCovered = mgr.constantFalse();
        BF strategy = mgr.constantFalse();
        for (auto it = strategyDumpingData.begin();
             it != strategyDumpingData.end(); it++) {
            if (it->first == i) {
                BF newCases = it->second.ExistAbstract(varCubePostOutput) &
                              !casesCovered;
                strategy |= newCases & it->second;
                casesCovered |= newCases;
            }
        }
        positionalStrategiesForTheIndividualGoals[i] = strategy;

        // std::ostringstream gsName;
        // gsName << "/tmp/generalStrategy" << i << ".dot";
        // BF dump = variables[4] & !variables[6]& !variables[8] & strategy;
        // BF_newDumpDot(*this,dump,"PreInput PreOutput PostInput
        // PostOutput",gsName.str().c_str());
    }

    // Allocate counter variables
    for (unsigned int i = 1; i <= livenessGuarantees.size(); i = i << 1) {
        std::ostringstream os;
        os << "_jx_b" << counterVarNumbersPre.size();
        counterVarNumbersPre.push_back(
            addVariable(SymbolicStrategyCounterVar, os.str()));
        if (systemGoalEncoded)
            counterVarNumbersPost.push_back(
                addVariable(SymbolicStrategyCounterVar, os.str() + "'"));
    }

    if (!systemGoalEncoded) {
        goalTransitionSelectorVar =
            addVariable(SymbolicStrategyCounterVar, "strat_type");
    }

    computeVariableInformation();

    BF combinedStrategy = mgr.constantFalse();
    for (unsigned int i = 0; i < livenessGuarantees.size(); i++) {
        BF thisEncoding = mgr.constantTrue();
        for (unsigned j = 0; j < counterVarNumbersPre.size(); j++) {
            if (i & (1 << j)) {
                thisEncoding &= variables[counterVarNumbersPre[j]];
            } else {
                thisEncoding &= !variables[counterVarNumbersPre[j]];
            }
        }

        if (systemGoalEncoded) {
            // SystemGoalEncoded -- The full format for deterministic
            // systems Here, we also include the transitions for moving to
            // the next goal.
            BF thisEncodingPostStay = mgr.constantTrue();
            for (unsigned j = 0; j < counterVarNumbersPre.size(); j++) {
                if (i & (1 << j)) {
                    thisEncodingPostStay &=
                        variables[counterVarNumbersPost[j]];
                } else {
                    thisEncodingPostStay &=
                        !variables[counterVarNumbersPost[j]];
                }
            }

            BF thisEncodingPostGo = mgr.constantTrue();
            for (unsigned j = 0; j < counterVarNumbersPre.size(); j++) {
                if (((i + 1) % livenessGuarantees.size()) & (1 << j)) {
                    thisEncodingPostGo &=
                        variables[counterVarNumbersPost[j]];
                } else {
                    thisEncodingPostGo &=
                        !variables[counterVarNumbersPost[j]];
                }
            }

            BF theseTrans =
                (thisEncodingPostGo & livenessGuarantees[i]) |
                (thisEncodingPostStay & !(livenessGuarantees[i]));
            combinedStrategy |=
                thisEncoding &
                positionalStrategiesForTheIndividualGoals[i] & theseTrans;

        } else {
            // not SystemGoalEncoded -- The LTLMoP format
            combinedStrategy |=
                thisEncoding &
                positionalStrategiesForTheIndividualGoals[i] &
                ((!variables[goalTransitionSelectorVar]) |
                 livenessGuarantees[i]);
        }
    }

    return positionalStrategiesForTheIndividualGoals;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::computeExplicitStrategy() {
    // synthesize repeatly until it passes the intermediate state check
    std::vector<std::vector<int>> regAsns;
    std::vector<BF> positionalStrategiesForTheIndividualGoals;
    std::vector<std::vector<int>> planPrefix;
    while (true) {
        // check unsafe set
        BF unsafeIntermediate = safetyTransitionSet();
        safetySys &= unsafeIntermediate;
        // Get symbolic strategy first
        T::execute();
        positionalStrategiesForTheIndividualGoals = computeSymbolicStrategy();
        // prefix: initial condition to the first goal
        BF prefixStrategy = positionalStrategiesForTheIndividualGoals[0];
        std::vector<int> initAsn = pcvz->getRegionInit();
        regAsns.clear();
        BF initState = regAsn2BF(initAsn);
        BF curState = initState;
        while (true) {
            std::vector<int> curRA;
            // discretize the state
            for (unsigned int i = 0; i < variables.size(); i++) {
                // (!) the sequence should be A, A' B, B' ...so the sequence will work
                if (doesVariableInheritType(i, Pre)) {
                    curRA.push_back((((curState & variables[i]).isFalse()) ? 0 : 1));
                }
            }
            regAsns.push_back(curRA);
            // check if the goal is reached
            if (!(curState & livenessGuarantees[0]).isFalse()) {
                break;
            }
            // calculate the next state
            curState &= positionalStrategiesForTheIndividualGoals[0];
            curState = determinize(curState, postVars);
            curState = curState.ExistAbstract(varCubePre).SwapVariables(varVectorPre, varVectorPost);
        }
        // Do the CP for that
        pcvz->printMiniZinc(toPatcher, planassignment, "prefixAsgn", -1, -1, regAsns, false);
        planPrefix = pcvz->getPatch();
        // check for transitions
        BF isFB = getIntermediateStateFeedback(planPrefix);   // Intermediate State Feed Back
        if (isFB.isTrue()) {
            break;
        }
        // std::cout << "Intermediate State Checker Failed. " << std::endl;
        safetySys &= isFB;
    }
    // add the plan to the patcher
    for (int i = 0; i < planPrefix.size(); i++) {
        p4Plan.addTransitionPrefix(planPrefix[i]);
    }
    // the last region assignment in the prefix
    std::vector<int> rA_firstGoal = pcvz->edge2RegionState(false, planPrefix[planPrefix.size() - 1]);
    std::vector<int> rA = rA_firstGoal;
    // suffix: first goal to second, second to third..., last to the first, and do the CP for each
    for (int g = 0; g < livenessGuarantees.size(); g++) {
        std::cout << "Generate Suffix Part of the plan" << std::endl;
        // update the new initial condition
        for (int i = 0; i < pcvz->getRegionNum(); i++) {
            pcvz->setRegionConstByIndex(i, -1, rA[i], -1);
        }
        // for the "last to the first", do also the final condition
        if (g == (livenessGuarantees.size() - 1)) {
            for (int i = 0; i < pcvz->getRegionNum(); i++) {
                pcvz->setRegionConstByIndex(i, -1, -1, rA_firstGoal[i]);
            }
        }
        // calculate the symbolic plan
        std::vector<std::vector<int>> regAsns_g;
        std::vector<std::vector<int>> planSuffix;
        while (true) {
            // synthesis explicit strategy
            int goalID = (g + 1) % livenessGuarantees.size();
            BF suffixStrategy_g = positionalStrategiesForTheIndividualGoals[goalID];
            std::vector<int> initAsn_g = pcvz->getRegionInit();
            regAsns_g.clear();
            BF initState_g = regAsn2BF(initAsn_g);
            BF curState_g = initState_g;
            while (true) {
                std::vector<int> curRA;
                // discretize the state
                for (unsigned int i = 0; i < variables.size(); i++) {
                    // (!) the sequence should be A, A' B, B' ...so the sequence will work
                    if (doesVariableInheritType(i, Pre)) {
                        curRA.push_back((((curState_g & variables[i]).isFalse()) ? 0 : 1));
                    }
                }
                regAsns_g.push_back(curRA);
                // check if the goal is reached
                // TODO: if the guarantee is a set, then the previous guarantee[0] don't have to be this one, require modification here
                if (!(curState_g & livenessGuarantees[goalID]).isFalse()) {
                    break;
                }
                // calculate the next state
                curState_g &= suffixStrategy_g;
                curState_g = determinize(curState_g, postVars);
                curState_g = curState_g.ExistAbstract(varCubePre).SwapVariables(varVectorPre, varVectorPost);
            }
            // Do the CP for that
            pcvz->printMiniZinc(toPatcher, planassignment, "prefixAsgn", -1, -1, regAsns_g, g == (livenessGuarantees.size() - 1));
            planSuffix = pcvz->getPatch();
            // check for intermediate states
            BF isFB = getIntermediateStateFeedback(planSuffix);   // Intermediate State Feed Back
            if (isFB.isTrue()) {
                break;
            }
            // std::cout << "Intermediate State Checker Failed. " << std::endl;
            safetySys &= isFB;
            // Get symbolic strategy
            T::execute();
            positionalStrategiesForTheIndividualGoals = computeSymbolicStrategy();
        }
        p4Plan.newGoalSuffix();
        for (int i = 0; i < planSuffix.size(); i++) {
            p4Plan.addTransitionSuffix(planSuffix[i], g);
        }
        rA = pcvz->edge2RegionState(false, planSuffix[planSuffix.size() - 1]);
    }
    // print a dot for testing
    std::vector<std::vector<int>> wholePlan = p4Plan.getAllTransitions();
    pcvz->printPatch2Dot("originalPlans", wholePlan);
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::reallocation() {
    std::vector<std::vector<int>> dum;
    if (reasgnRegionState.size() != 0) {
        std::cout << "=========<Do the Re-Assignment: Should not work with Changing the capacity or Removing the edge>==========" << std::endl;
        // Check if the reassignment value is valid
        std::vector<int> maxNs = pcvz->getRegionCap();
        auto tgtTransAs = p4Plan.getTransition(reasgnGoalNum, reasgnStateNum);
        std::vector<int> tgtStateAs;
        // check if we want the last state in the pre/suffix
        if (tgtTransAs.size() == 0) {
            tgtTransAs = p4Plan.getTransition(reasgnGoalNum, reasgnStateNum - 1);
            tgtStateAs = pcvz->edge2RegionState(false, tgtTransAs);
        } else {
            tgtStateAs = pcvz->edge2RegionState(true, tgtTransAs);
        }
        for (int i = 0; i < maxNs.size(); i++) {
            if (reasgnRegionState[i] > maxNs[i]) {
                std::cout << "Re-assignment " << reasgnRegionState[i] << " is larger than the capacity" << maxNs[i] << std::endl;
                return;
            }
            if ((reasgnRegionState[i] == 0) != (tgtStateAs[i] == 0)) {
                std::cout << "Re-assignment " << reasgnRegionState[i] << " and original assignment " << tgtStateAs[i] << " have different state." << std::endl;
                return;
            }
        }

        // iteratively check for intermediate states
        std::vector<std::vector<int>> ppatch;
        BF isFB = safetyTransitionSet();
        while (true) {
            // calculate the fixpoint for the state
            BF raTranRM = getRMSG();
            if (isFB.isTrue()) {
                std::cout << "isFB is true" << std::endl;
            }
            BF rSG = safetySysNoRM & raTranRM & isFB;
            BFFixedPoint rMu1(mgr.constantFalse());
            bool checkReach = false;
            BF raState = regAsn2BF(tgtStateAs);
            BF lastLayer = raState;
            BF rstrategy = mgr.constantFalse();
            int layerNumb = 0;
            while (!(rMu1.isFixedPointReached())) {
                layerNumb++;
                BF pathFnd = rSG & (lastLayer.SwapVariables(varVectorPre, varVectorPost));
                BF thisLayer = pathFnd.ExistAbstract(varCubePostOutput);
                lastLayer |= thisLayer;
                // rstrategy |= pathFnd | (thisLayer & thisLayer.SwapVariables(varVectorPre, varVectorPost));
                rstrategy |= pathFnd;
                rMu1.update(lastLayer);
            }
            std::cout << "Max Re-assignment Patch Layer Number: " << layerNumb << std::endl;

            // do the CP, expand the patch size until the max level is reached
            int tranNumb = 1;
            for (int i = 0; i < maxNs.size(); i++) {
                pcvz->setRegionConstByIndex(i, -1, tgtStateAs[i], reasgnRegionState[i]);
            }
            pcvz->clearClauses();
            BF2CVZClauses(rstrategy, firstStrategy);
            while (tranNumb <= layerNumb) {
                pcvz->printMiniZinc(toPatcher, patching, "swarmTestReasgnFork", tranNumb, -1, dum, true);
                ppatch = pcvz->getPatch();
                if (ppatch.size() > 0) {
                    // pcvz->printPatch2Dot("ReAssignment");
                    std::cout << "patch Found" << std::endl;
                    break;
                }
                tranNumb++;
            }
            if (ppatch.size() == 0) {
                std::cout << "Weird... Reallocation can't found." << std::endl;
                return;
            }
            // check for intermediate states
            BF isFB_this = getIntermediateStateFeedback(ppatch);   // Intermediate State Feed Back
            isFB &= isFB_this;
            if (isFB_this.isTrue()) {
                break;
            }
        }

        // update the new assignment to the suffix part of the original patch
        // will have to use another CP...
        // new prefix
        std::vector<std::vector<int>> oldAsgn;
        for (int i = reasgnStateNum; i < p4Plan.getSuffixSize(reasgnGoalNum); i++) {
            oldAsgn.push_back(p4Plan.getTransition(reasgnGoalNum, i));
        }
        std::cout << "..." << std::endl;
        std::vector<int> rA_firstGoal;
        std::vector<std::vector<int>> prefixPatch;
        if (oldAsgn.size() != 0) {
            for (int i = 0; i < maxNs.size(); i++) {
                // pcvz->setRegionConstByIndex(i, -1, reasgnRegionState[i], plan[plan.size() - 1][i]);
                pcvz->setRegionConstByIndex(i, -1, reasgnRegionState[i], -1);
            }
            pcvz->printMiniZinc(toPatcher, reassignment, "newAssignPrefix", -1, -1, oldAsgn, false);
            prefixPatch = pcvz->getPatch();
            if (prefixPatch.size() == 0) {
                std::cout << "Update the plan failed after the re-allocation. " << std::endl;
                return;
            }
            rA_firstGoal = pcvz->edge2RegionState(false, prefixPatch.back());
        } else {
            rA_firstGoal = pcvz->edge2RegionState(false, ppatch.back());
        }
        std::vector<int> rA = rA_firstGoal;
        std::vector<std::vector<std::vector<int>>> suffixes;
        // update the plan with the new prefix
        // new suffix
        for (int g = 0; g < livenessGuarantees.size(); g++) {
            int goalID = (g + reasgnGoalNum + 1) % livenessGuarantees.size();
            oldAsgn.clear();
            for (int i = 0; i < p4Plan.getSuffixSize(goalID); i++) {
                oldAsgn.push_back(p4Plan.getTransition(goalID, i));
            }
            for (int i = 0; i < maxNs.size(); i++) {
                // rA_firstGoal is only used when g == reasgnGoalNum
                pcvz->setRegionConstByIndex(i, -1, rA[i], rA_firstGoal[i]);
            }
            pcvz->printMiniZinc(toPatcher, reassignment, "newAssignSuffix", -1, -1, oldAsgn, (g == (livenessGuarantees.size() - 1)));
            auto suffixPatch = pcvz->getPatch();
            if (suffixPatch.size() == 0) {
                std::cout << "Update the plan failed after the re-allocation. " << std::endl;
                return;
                // break;
            }
            rA = pcvz->edge2RegionState(false, suffixPatch.back());
            suffixes.push_back(suffixPatch);
        }
        // make patch if the updated assignment exceed capacities somewhere else
        p4Plan.clean();
        for (int i = 0; i < ppatch.size(); i++) {
            p4Plan.addTransitionPrefix(ppatch[i]);
        }
        if (oldAsgn.size() != 0) {
            for (int i = 0; i < prefixPatch.size(); i++) {
                p4Plan.addTransitionPrefix(prefixPatch[i]);
            }
        }
        for (int g = 0; g < suffixes.size(); g++) {
            p4Plan.newGoalSuffix();
            for (int i = 0; i < suffixes[g].size(); i++) {
                p4Plan.addTransitionSuffix(suffixes[g][i], g);
            }
        }

        // print as dot
        pcvz->printPatch2Dot("ReAssigned", p4Plan.getAllTransitions());

        // update p4Plan edge capacity
        std::vector<int> edgeMax = pcvz->getEdgeMax();
        std::vector<int> edgeIDs;
        for (int i = 0; i < edgeMax.size(); i++) {
            edgeIDs.push_back(i);
        }
        p4Plan.updateCapacity(edgeIDs, edgeMax);
        // return;
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::updateRAS() {
    rAs.clear();
    // p4Plan.printTransitions();
    for (int i = 0; i <= livenessGuarantees.size(); i++) {
        std::vector<std::vector<int>> trans;
        for (int j = 0; j < p4Plan.getSuffixSize(i - 1); j++) {
            trans.push_back(pcvz->edge2RegionState(true, p4Plan.getTransition(i - 1, j)));
        }
        trans.push_back(pcvz->edge2RegionState(false, p4Plan.getTransition(i - 1, p4Plan.getSuffixSize(i - 1) - 1)));
        rAs.push_back(trans);
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
std::pair<BF, int>
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::checkReachable(BF start, BF end, BF &isFB) {
    // transition rules of the region map
    BF newTranRM = getRMSG();
    // the overall safety guarantees
    BF pSG = safetySysNoRM & newTranRM & isFB;
    // initialize the checker
    BFFixedPoint pMu1(mgr.constantFalse());
    bool checkReach = false;
    BF lastLayer = end;
    BF pstrategy = mgr.constantFalse();
    int layerNumb = 0;
    // fix point operation starts
    while (true) {
        if (pMu1.isFixedPointReached()) {
            break;
        }
        layerNumb++;
        BF pathFnd = pSG & (lastLayer.SwapVariables(varVectorPre, varVectorPost));
        BF thisLayer = pathFnd.ExistAbstract(varCubePostOutput);
        if (!(thisLayer & start).isFalse()) {
            checkReach = true;
        }
        lastLayer |= thisLayer;
        // pstrategy |= pathFnd | (thisLayer & thisLayer.SwapVariables(varVectorPre, varVectorPost));
        pstrategy |= pathFnd;
        pMu1.update(lastLayer);
    }
    return std::make_pair(pstrategy, layerNumb);
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
std::vector<std::vector<int>>
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::patchForGivenHorizon(std::vector<int> iniStateAs, std::vector<int> fnlStateAs, int expLayNumb, bool final) {
    BF iniState = regAsn2BF(iniStateAs);
    BF fnlState = regAsn2BF(fnlStateAs);
    // iteratively find a patch that does not violate the intermediate states
    BF isFB = safetyTransitionSet();   // mgr.constantTrue();
    std::vector<std::vector<int>> ppatch;
    while (true) {
        // check reachability: one goal, one init and no env. vars...
        std::cout << "=========<Patching Iteration: Check Reachability>==========" << std::endl;
        std::pair<BF, int> reachablePair = checkReachable(iniState, fnlState, isFB);
        BF pstrategy = reachablePair.first;
        int layerNumb = reachablePair.second;
        bool checkReach = (layerNumb != 0);
        std::cout << "Patch Layer Number Guess: " << layerNumb << std::endl;
        // if reachable, synthesis the strategy and do the cp
        if (checkReach) {
            std::cout << "=========<Patching Iteration: Trying CP the Patch>==========" << std::endl;
            for (int i = 0; i < iniStateAs.size(); i++) {
                pcvz->setRegionConstByIndex(i, -1, iniStateAs[i], fnlStateAs[i]);
            }

            pcvz->clearClauses();
            BF2CVZClauses(pstrategy, firstStrategy);
        }
        // increase the layer number until a solution is found
        ppatch.clear();
        for (int i = expLayNumb; i <= layerNumb; i++) {
            std::cout << "Trying for " << i << " layers" << std::endl;
            std::vector<std::vector<int>> dum;
            pcvz->printMiniZinc(toPatcher, patching, "swarmTestFork", i, -1, dum, final);
            ppatch = pcvz->getPatch();
            if (ppatch.size() > 0) {
                break;
            }
        }
        if (ppatch.size() > 0) {
            // check for intermediate states
            BF isFB_this = getIntermediateStateFeedback(ppatch);   // Intermediate State Feed Back
            isFB &= isFB_this;
            if (isFB_this.isTrue()) {
                std::cout << "patch Found" << std::endl;
                return ppatch;
            }
        } else {
            // synthesis failed for the given initial/final state
            return ppatch;
        }
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
std::pair<std::vector<std::vector<int>>, int>
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::patchForGivenHorizonDoubleStrategy(std::vector<int> iniStateAs, std::vector<int> fnlStateAs, int expFirstLayNumb, int expSecondLayNumb, int goalNumb, bool final) {
    std::vector<std::vector<int>> ppatch;
    auto ret = std::make_pair(ppatch, 0);
    // specify the goal state within the patch
    BF goal;
    if ((goalNumb < livenessGuarantees.size()) && (goalNumb >= 0)) {
        goal = livenessGuarantees[goalNumb];
    } else {
        std::cout << "Liveness Guarantee Index Out of Range: " << goalNumb << std::endl;
        return ret;
    }
    // initialize variables
    BF iniState = regAsn2BF(iniStateAs);
    BF fnlState = regAsn2BF(fnlStateAs);
    // iteratively find a patch that does not violate the intermediate states
    BF isFB = safetyTransitionSet();   // mgr.constantTrue();
    while (true) {
        // check reachability for initStarte to the goal: one goal, one init and no env. vars...
        std::cout << "=========<Patching Iteration: Check Reachability>==========" << std::endl;
        std::pair<BF, int> reachablePair = checkReachable(iniState, goal, isFB);
        BF strategy2Goal = reachablePair.first;
        int layerNumb2Goal = reachablePair.second;
        bool checkReach2Goal = (layerNumb2Goal != 0);
        // check reachability from goal to the final state: one goal, one init and no env. vars...
        reachablePair = checkReachable(goal, fnlState, isFB);
        BF strategy2Fnl = reachablePair.first;
        int layerNumb2Fnl = reachablePair.second;
        bool checkReach2Fnl = (layerNumb2Fnl != 0);
        if ((!checkReach2Goal) || (!checkReach2Fnl)) {
            std::cout << "This horizon is not reachable" << std::endl;
            return ret;
        }
        // if reachable, synthesis the strategy and do the cp
        std::cout << "Patch Layer Number Guess: " << layerNumb2Goal + layerNumb2Fnl << std::endl;
        std::cout << "=========<Patching Iteration: Trying CP the Patch>==========" << std::endl;
        for (int i = 0; i < iniStateAs.size(); i++) {
            pcvz->setRegionConstByIndex(i, -1, iniStateAs[i], fnlStateAs[i]);
        }
        pcvz->clearClauses();
        BF2CVZClauses(strategy2Goal, firstStrategy);
        BF2CVZClauses(strategy2Fnl, secondStrategy);
        BF2CVZClauses(goal, goalState);
        // increase the layer number until a solution is found
        ppatch.clear();
        int maxLayerNumbAdd = ((layerNumb2Goal - expFirstLayNumb) > (layerNumb2Fnl - expSecondLayNumb)) ? (layerNumb2Goal - expFirstLayNumb) : (layerNumb2Fnl - expSecondLayNumb);
        int i = 0;
        for (i = 0; i <= maxLayerNumbAdd; i++) {
            std::cout << "Trying for " << i << " layers" << std::endl;
            std::vector<std::vector<int>> dum;
            pcvz->printMiniZinc(toPatcher, doublestrategy, "swarmTestFork", expFirstLayNumb + i, expSecondLayNumb + i, dum, final);
            ppatch = pcvz->getPatch();
            if (ppatch.size() > 0) {
                break;
            }
        }
        // return the patch
        if (ppatch.size() > 0) {
            // check for intermediate states
            BF isFB_this = getIntermediateStateFeedback(ppatch);   // Intermediate State Feed Back
            isFB &= isFB_this;
            if (isFB_this.isTrue()) {
                std::cout << "patch Found" << std::endl;
                ret.first = ppatch;
                ret.second = expFirstLayNumb + i;
                return ret;
            }
        } else {
            // synthesis failed for the given initial/final state
            return ret;
        }
    }
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
bool
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::patchForGoal_Legacy(int goalID, std::vector<std::pair<int, int>> locP) {
    std::cout << "=========<Patching For Goal: " << goalID << " >==========" << std::endl;
    std::vector<std::vector<int>> dum;
    // check for the patch bound
    int tMax = p4Plan.getSuffixSize(goalID) - 1;
    int tIni = tMax, tFnl = 0;
    for (int i = 0; i < locP.size(); i++) {
        if (locP[i].first == goalID) {
            int tN = locP[i].second;
            std::cout << "Problematic State: " << tN << std::endl;
            if (tN > tFnl) {
                tFnl = tN;
            }
            if (tN < tIni) {
                tIni = tN;
            }
        }
    }
    // do a while loop until we find a patch
    std::cout << "=========<Patching Iteration>==========" << std::endl;
    bool patchFound = false;
    if (tIni > tFnl) {
        std::cout << "No need to modify in this case" << std::endl;
        patchFound = true;
        return patchFound;
    }
    BF isFB = mgr.constantTrue();
    std::vector<std::vector<int>> ppatch;
    while (!patchFound) {
        // get init/goal state as BF
        std::cout << "=========<Patching Iteration: Getting Boundaries>==========" << std::endl;
        std::cout << "Initial State: " << tIni << ", final state: " << tFnl << std::endl;
        // auto iniStateAs = pcvz->edge2RegionState(true, p4Plan.getTransition(goalID, tIni));
        // auto fnlStateAs = pcvz->edge2RegionState(false, p4Plan.getTransition(goalID, tFnl));
        auto iniStateAs = rAs[goalID + 1][tIni];
        auto fnlStateAs = rAs[goalID + 1][tFnl + 1];
        BF iniState = regAsn2BF(iniStateAs);
        BF fnlState = regAsn2BF(fnlStateAs);
        // Build the patch
        // Note that if the final state is the goal state, the final numeric condition is not important
        ppatch = patchForGivenHorizon(iniStateAs, fnlStateAs, (tFnl - tIni + 1), (tFnl != tMax));
        if (ppatch.size() > 0) {
            p4Plan.makePatch(ppatch, goalID, tIni, tFnl);
            patchFound = true;
            break;
        }
        // if not reachable or the cp failed, increase the patch horizon
        std::cout << "=========<Patching Iteration: Expand the Patch Horizon>==========" << std::endl;
        int ntIni = std::max(0, tIni - 1);
        int ntFnl = std::min((int) tMax, tFnl + 1);
        if ((ntIni == tIni) && (ntFnl == tFnl)) {
            break;
        } else {
            tIni = ntIni;
            tFnl = ntFnl;
        }
    }

    if (!patchFound) {
        return patchFound;
    }

    // if the goal state is modified, add a patch to the next transition
    if (tFnl == tMax) {
        std::cout << "=========<Patching Extension: Goal State Modified>==========" << std::endl;
        std::vector<int> goalStateAs = pcvz->edge2RegionState(false, ppatch.back());
        // iteratively find the patch horizon
        int nextGoalID = (goalID + 1) % livenessGuarantees.size();
        int nextSuffixSize = p4Plan.getSuffixSize(nextGoalID);   // the # of transitions in the next suffix
        // we don't consider the patch that covers the next goal state
        for (int i = 0; i < nextSuffixSize; i++) {
            // try to patch
            std::vector<int> fnlStateAs = rAs[nextGoalID + 1][i + 1];
            std::vector<std::vector<int>> spatch = patchForGivenHorizon(goalStateAs, fnlStateAs, i + 1, true);
            if (spatch.size() > 0) {
                p4Plan.makePatch(spatch, nextGoalID, 0, i);
                // And also patch for the next rAs
                std::vector<std::vector<int>> rAs_ng_old = rAs[nextGoalID + 1];
                rAs[nextGoalID + 1].clear();
                for (int j = 0; j < spatch.size(); j++) {
                    rAs[nextGoalID + 1].push_back(pcvz->edge2RegionState(true, spatch[j]));
                }
                for (int j = i + 1; j < rAs_ng_old.size(); j++) {
                    rAs[nextGoalID + 1].push_back(rAs_ng_old[j]);
                }
                break;
            }
            // if not found until i = nSS - 1, then the patching is failed
            if (i == nextSuffixSize - 1) {
                patchFound = false;
            }
        }
    }

    return patchFound;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
bool
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::patchForGoal(int goalID, std::vector<std::pair<int, int>> locP) {
    std::cout << "=========<Patching For Goal: " << goalID << " >==========" << std::endl;
    std::vector<std::vector<int>> dum;
    // check for the patch bound
    int nextGoalID = (goalID + 1) % livenessGuarantees.size();
    int tMax = p4Plan.getSuffixSize(goalID) - 1;
    int tMaxNext = p4Plan.getSuffixSize(nextGoalID) - 1;
    int tIni = tMax, tFnl = 0;
    for (int i = 0; i < locP.size(); i++) {
        if (locP[i].first == goalID) {
            int tN = locP[i].second;
            std::cout << "Problematic State: " << tN << std::endl;
            if (tN > tFnl) {
                tFnl = tN;
            }
            if (tN < tIni) {
                tIni = tN;
            }
        }
    }
    // do a while loop until we find a patch
    bool doubleStrategy = false;
    std::cout << "=========<Patching Iteration>==========" << std::endl;
    bool patchFound = false;
    if (tIni > tFnl) {
        std::cout << "No need to modify in this case" << std::endl;
        patchFound = true;
        return patchFound;
    }
    BF isFB = safetyTransitionSet();   // mgr.constantTrue();
    std::vector<std::vector<int>> ppatch;
    while (!patchFound) {
        // get init/goal state as BF
        std::cout << "=========<Patching Iteration: Getting Boundaries>==========" << std::endl;
        std::cout << "Initial State: " << tIni << ", final state: " << tFnl << std::endl;
        auto iniStateAs = rAs[goalID + 1][tIni];
        // in the double strategy case, the final state comes from the next suffix
        std::vector<int> fnlStateAs;
        if (doubleStrategy) {
            fnlStateAs = rAs[nextGoalID + 1][tFnl + 1];
        } else {
            fnlStateAs = rAs[goalID + 1][tFnl + 1];
        }
        BF iniState = regAsn2BF(iniStateAs);
        BF fnlState = regAsn2BF(fnlStateAs);
        // Build the patch
        // Note that if the final state is the goal state, the final numeric condition is not important
        if (doubleStrategy) {
            int goalTranID = tMax;
            // pair := (patch, the_#trans_of_first_strategy)
            auto ppatchPair = patchForGivenHorizonDoubleStrategy(iniStateAs, fnlStateAs, (tMax - tIni + 1), (tFnl + 1), (goalID + 1) % livenessGuarantees.size(), true);
            ppatch = ppatchPair.first;
            if (ppatch.size() > 0) {
                int fsLen = ppatchPair.second;   // first strategy length
                int ssLen = ppatch.size() - ppatchPair.second;
                std::vector<std::vector<int>> firstPatch(ppatch.begin(), ppatch.begin() + fsLen);
                std::vector<std::vector<int>> secondPatch(ppatch.begin() + fsLen, ppatch.begin() + ppatch.size());
                p4Plan.makePatch(firstPatch, goalID, tIni, tMax);
                p4Plan.makePatch(secondPatch, nextGoalID, 0, tFnl);
                // And also patch for the next rAs
                std::vector<std::vector<int>> rAs_ng_old = rAs[nextGoalID + 1];
                rAs[nextGoalID + 1].clear();
                for (int j = 0; j < ssLen; j++) {
                    rAs[nextGoalID + 1].push_back(pcvz->edge2RegionState(true, ppatch[fsLen + j]));
                }
                for (int j = tFnl + 1; j < rAs_ng_old.size(); j++) {
                    rAs[nextGoalID + 1].push_back(rAs_ng_old[j]);
                }
                patchFound = true;
                break;
            }
        } else {
            ppatch = patchForGivenHorizon(iniStateAs, fnlStateAs, (tFnl - tIni + 1), true);
            if (ppatch.size() > 0) {
                p4Plan.makePatch(ppatch, goalID, tIni, tFnl);
                patchFound = true;
                break;
            }
        }
        // if not reachable or the cp failed, increase the patch horizon
        std::cout << "=========<Patching Iteration: Expand the Patch Horizon>==========" << std::endl;
        int ntIni = std::max(0, tIni - 1);
        int ntFnl;
        if (doubleStrategy) {
            ntFnl = std::min((int) tMaxNext, tFnl + 1);
        } else {
            ntFnl = std::min((int) tMax, tFnl + 1);
        }
        // this means the patch will cover the goal state
        if ((!doubleStrategy) && (ntFnl == tFnl)) {
            doubleStrategy = true;
            ntFnl = 0;
        }
        // max patch horizon will be "the start of the current suffix" and "the end of the next suffix"
        if ((ntIni == tIni)) {
            if ((doubleStrategy) && (ntFnl == tFnl)) {
                break;
            }
        }
        tIni = ntIni;
        tFnl = ntFnl;
    }

    return patchFound;
}

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::computeAndPrintSymbolicStrategy(std::string filename) {

    // We don't want any reordering from this point onwards, as
    // the BDD manipulations from this point onwards are 'kind of simple'.
    mgr.setAutomaticOptimisation(false);

    // before synthesis, check the layer
    // updateL2G();

    // safetySys2IntState();
    std::cout << "=========<Generate an original plan by CP>==========" << std::endl;
    computeExplicitStrategy();
    std::vector<std::vector<int>> dum;
    // updateL2G();

    while (true) {
        std::cout << "=========<Start building the patch>==========" << std::endl;

        // std::string inp;
        // std::cout << "Type something..." << std::endl;
        // std::cin >> inp;
        // std::cout << inp << std::endl;

        // input for specify modification
        std::string typeM, modMaxStr, gidStr, tidStr, asnStr;
        rmvEdgeFrom = rmvEdgeTo = modReg = "";
        reasgnRegionState.clear();
        std::cout << "Please specify the modification type: {DEC/RMV/RAL}" << std::endl;
        std::cin >> typeM;
        if (typeM == "DEC") {
            std::cout << "The region name?" << std::endl;
            std::cin >> modReg;
            std::cout << "And the new capacity?" << std::endl;
            std::cin >> modMaxStr;
            modMax = stoi(modMaxStr);
        } else if (typeM == "RMV") {
            std::cout << "Remove the edge.  From which region?" << std::endl;
            std::cin >> rmvEdgeFrom;
            std::cout << "To which region?" << std::endl;
            std::cin >> rmvEdgeTo;
        } else if (typeM == "RAL") {
            std::cout << "Please specify the transition ID: " << std::endl;
            std::cin >> tidStr;
            std::cout << "and the goal ID: " << std::endl;
            std::cin >> gidStr;
            reasgnStateNum = stoi(tidStr);
            reasgnGoalNum = stoi(gidStr);
            std::cout << "and the expected number assignment on each region: " << std::endl;
            for (int i = 0; i < pcvz->getRegionNum(); i++) {
                std::cout << pcvz->getRegionName(i) << ": " << std::endl;
                std::cin >> asnStr;
                reasgnRegionState.push_back(stoi(asnStr));
            }
        } else {
            std::cout << "Strange Modification Type " << typeM << std::endl;
        }

        // check for reassignment
        reallocation();
        updateRAS();

        std::cout << "=========<Do an initial guess on Patch Horizon>==========" << std::endl;
        // If there is any removal of an edge
        if ((rmvEdgeFrom.size() != 0) && (rmvEdgeTo.size() != 0)) {
            bool directed = false;
            if (rmvEdgeFrom == rmvEdgeTo) {
                directed = true;
                pcvz->addEdgeByName(directed, true, rmvEdgeFrom, rmvEdgeTo);
                int rmvEdgeID = pcvz->getEdgeIDByName(rmvEdgeFrom, rmvEdgeTo);
                p4Plan.updateEdge(true, rmvEdgeID);
            } else {
                directed = true;
                int rmvEdgeID1 = pcvz->getEdgeIDByName(rmvEdgeFrom, rmvEdgeTo);
                pcvz->addEdgeByName(directed, true, rmvEdgeFrom, rmvEdgeTo);
                p4Plan.updateEdge(true, rmvEdgeID1);
                int rmvEdgeID2 = pcvz->getEdgeIDByName(rmvEdgeTo, rmvEdgeFrom);
                pcvz->addEdgeByName(directed, true, rmvEdgeTo, rmvEdgeFrom);
                p4Plan.updateEdge(true, rmvEdgeID2);
                // p4Plan.cleanEidRmv();
                // p4Plan.specifyEidRmv(rmvEdgeID1);
                // p4Plan.specifyEidRmv(rmvEdgeID2);
            }
            // update the max region capaciy
            auto newEdgeMax = pcvz->getEdgeMax();
            std::vector<int> edgeIDs;
            for (int i = 0; i < newEdgeMax.size(); i++) {
                std::cout << newEdgeMax[i] << " ";
                edgeIDs.push_back(i);
            }
            std::cout << std::endl;
            p4Plan.updateCapacity(edgeIDs, newEdgeMax);
            p4Plan.printTransitions();
        }
        // If the modification of regional robot capacity is required
        std::vector<int> originalEdgeMax, newEdgeMax;
        // std::cout << modReg << std::endl;
        if ((modReg.size() != 0)) {
            originalEdgeMax = pcvz->getEdgeMax();
            pcvz->setRegionConstByName(modReg, modMax, -1, -1);
            newEdgeMax = pcvz->getEdgeMax();
            // update the max region capaciy
            std::vector<int> edgeIDs;
            for (int i = 0; i < newEdgeMax.size(); i++) {
                std::cout << newEdgeMax[i] << " ";
                edgeIDs.push_back(i);
            }
            std::cout << std::endl;
            p4Plan.updateCapacity(edgeIDs, newEdgeMax);
        }

        std::vector<std::pair<int, int>> locP = p4Plan.localizePatch();

        // do the patch
        bool patchFound = true;
        for (int i = 0; i <= livenessGuarantees.size(); i++) {
            int goalID = i - 1;
            locP = p4Plan.localizePatch();
            bool isFound = patchForGoal(goalID, locP);
            patchFound &= isFound;
            if (isFound) {
                std::cout << "Goal " << i << " is patched" << std::endl;
            } else {
                std::cout << "Goal " << i << " patch fails" << std::endl;
            }
        }

        p4Plan.printTransitions();
        // if failed, print something
        if (!patchFound) {
            std::cout << "Patch Failed" << std::endl;
        } else {
            // print as dot
            pcvz->printPatch2Dot("Patched", p4Plan.getAllTransitions());
        }
    }

    // TODO Lists:
    // 1. Integrate Ji's and Salar's Plan Synthesizer (\/)
    // 2. Make current codes be able to run plan with many goals (\/)
    // 3. Integrate Salar's Intermediate State Checker (\/)
    // 3-1. Patching prefix after reallocation causes problems because the patch is difficult to include both the re-allocation and new prefixes
    // 4. DO the patch when the modification affects the goal states (\/)
    // 5. Interactive interface (\/)
    // 6. Difference of the BB in CP solvers and doing "feedback and conflict solving" on GR(1)?
    // 7. Cases:
    // 7-1. DEC O 6 -> RMV M R: Runs forever
    // 7-2. RMV P O -> DEC P 6: Runs OK
    // 7-3. RAL 1 -1 (4,6) ->6: Runs OK
    // 8. Patch Expansion
    // 9. Cost Function When Reallocate (\/)
    // 10. the "goal state" should be the same for prefix and last suffix (both go to live. guar. 0) after any modification
}

template class XSwarmTest<GR1Context, false, false>;
