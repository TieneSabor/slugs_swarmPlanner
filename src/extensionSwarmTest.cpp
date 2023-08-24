#include <extensionSwarmTest.hpp>

using namespace swarmTest;

template <class T, bool oneStepRecovery, bool systemGoalEncoded>
void
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::init(std::list<std::string> &filenames) {
    // build the regional map
    BF RMSG = parseRegionMap(filenames);
    // std::cout << "?" << std::endl;
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
    // std::cout << variables.size() << std::endl;
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
        } else if (0 == (currentLine.compare("# regional map end"))) {
            readMode = 0;
            break;
        } else {
            std::istringstream is(currentLine);
            std::string shrp, reg, cap, ini, fnl, from, to;
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
XSwarmTest<T, oneStepRecovery, systemGoalEncoded>::computeAndPrintSymbolicStrategy(std::string filename) {

    // We don't want any reordering from this point onwards, as
    // the BDD manipulations from this point onwards are 'kind of simple'.
    mgr.setAutomaticOptimisation(false);

    // before synthesis, check the layer
    // updateL2G();

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

    // make strategy BDD a set of CNF clauses
    int goalNumb = 0;
    combinedStrategy = positionalStrategiesForTheIndividualGoals[goalNumb];
    filename = filename + "_cnf";
    std::string header_cnf = "";
    int clauseBufferSize = 100000;
    int clausesInCNFIndex[clauseBufferSize], clauseN = 0, varSize = 0;
    mgr.writeBDDToCNFFile(filename.c_str(), header_cnf,
                          combinedStrategy, variables, variableNames, clausesInCNFIndex, clauseN, varSize);

    // Build miniZinc Converter
    std::cout << "=========<Add strategy as clauses>==========" << std::endl;
    // add clauses to the converter
    for (int i = 0; i < clauseN; i++) {
        pcvz->newClause();
        // std::cout << "Clause #" << i << " ";
        for (int j = 0; j < varSize; j++) {
            // std::cout << clausesInCNFIndex[i * varSize + j] << " ";
            std::string lit = variableNames[j];
            switch (clausesInCNFIndex[i * varSize + j]) {
            case 0:
                pcvz->addLiteral2LastClauseByName(lit);
                break;
            case 1:
                pcvz->addLiteral2LastClauseByName("-" + lit);
                break;
            case 2:
                break;
            default:
                std::cout << "Unknown literal value: " << clausesInCNFIndex[i * varSize + j] << std::endl;
            }
        }
        // std::cout << std::endl;
    }
    // print to file or terminal
    // pcvz->printMiniZinc(toTerminal, "swarmTest", 2);
    // solve the CP for only one liveness guarantee
    std::cout << "=========<Generate an original plan by CP>==========" << std::endl;
    updateL2G();
    int ly = checkLayer(initState, 0);
    std::cout << "Init state is at " << ly << "th layer" << std::endl;
    pcvz->printMiniZinc(toPatcher, "swarmTestFork", ly);
    pcvz->printPatch2Dot("OriginalPlan");
    auto plan = pcvz->getPatch();
    // build patcher for the one live. guar. case
    std::cout << "=========<Start building the patch>==========" << std::endl;
    patcher p4Plan;
    p4Plan.newGoalSuffix();
    for (int i = 0; i < plan.size(); i++) {
        p4Plan.addTransitionSuffix(plan[i], 0);
    }
    // p4Plan.printTransitions();
    // locate the patch start-end
    std::cout << "=========<Do an initial guess on Patch Horizon>==========" << std::endl;
    auto originalEdgeMax = pcvz->getEdgeMax();
    pcvz->setRegionConstByName(modReg, modMax, -1, -1);
    auto newEdgeMax = pcvz->getEdgeMax();
    std::vector<int> edgeIDs;
    for (int i = 0; i < newEdgeMax.size(); i++) {
        edgeIDs.push_back(i);
    }
    auto locP = p4Plan.localizePatchModMax(edgeIDs, newEdgeMax);
    int tIni = plan.size() - 1, tFnl = 0;
    for (int i = 0; i < locP.size(); i++) {
        int tN = locP[i].second;
        if (tN > tFnl) {
            tFnl = tN;
        }
        if (tN < tIni) {
            tIni = tN;
        }
    }
    // do a while loop until we find a patch
    std::cout << "=========<Patching Iteration>==========" << std::endl;
    bool patchFound = false;
    if (tIni > tFnl) {
        std::cout << "No need to modify (?)" << std::endl;
        patchFound = true;
    }
    while (!patchFound) {
        // get init/goal state as BF
        std::cout << "=========<Patching Iteration: Getting Boundaries>==========" << std::endl;
        auto iniStateAs = pcvz->edge2RegionState(true, p4Plan.getTransition(0, tIni));
        auto fnlStateAs = pcvz->edge2RegionState(false, p4Plan.getTransition(0, tFnl));
        BF iniState = regAsn2BF(iniStateAs);
        BF fnlState = regAsn2BF(fnlStateAs);
        // check reachability: one goal, one init and no env. vars...
        std::cout << "=========<Patching Iteration: Check Reachability>==========" << std::endl;
        BF newTranRM = getRMSG();
        BF pSG = safetySysNoRM & newTranRM;
        BFFixedPoint pMu1(mgr.constantFalse());
        bool checkReach = false;
        BF lastLayer = fnlState;
        BF pstrategy = mgr.constantFalse();
        int layerNumb = 0;
        while (!checkReach) {
            if (pMu1.isFixedPointReached()) {
                break;
            }
            layerNumb++;
            BF pathFnd = pSG & (lastLayer.SwapVariables(varVectorPre, varVectorPost));
            BF thisLayer = pathFnd.ExistAbstract(varCubePostOutput);
            if (!(thisLayer & iniState).isFalse()) {
                checkReach = true;
            }
            lastLayer |= thisLayer;
            pstrategy |= pathFnd;
            pMu1.update(lastLayer);
        }
        // if reachable, synthesis the strategy and do the cp
        if (checkReach) {
            std::cout << "=========<Patching Iteration: Trying CP the Patch>==========" << std::endl;
            for (int i = 0; i < iniStateAs.size(); i++) {
                pcvz->setRegionConstByIndex(i, -1, iniStateAs[i], fnlStateAs[i]);
            }
            filename = filename + "_cnf";
            std::string header_cnf = "";
            int pclausesInCNFIndex[1000000], pclauseN = 0, pvarSize = 0;
            mgr.writeBDDToCNFFile(filename.c_str(), header_cnf,
                                  pstrategy, variables, variableNames, pclausesInCNFIndex, pclauseN, pvarSize);
            pcvz->clearClauses();
            for (int i = 0; i < clauseN; i++) {
                pcvz->newClause();
                // std::cout << "Clause #" << i << " ";
                for (int j = 0; j < varSize; j++) {
                    // std::cout << clausesInCNFIndex[i * varSize + j] << " ";
                    std::string lit = variableNames[j];
                    switch (clausesInCNFIndex[i * varSize + j]) {
                    case 0:
                        pcvz->addLiteral2LastClauseByName(lit);
                        break;
                    case 1:
                        pcvz->addLiteral2LastClauseByName("-" + lit);
                        break;
                    case 2:
                        break;
                    default:
                        std::cout << "Unknown literal value: " << clausesInCNFIndex[i * varSize + j] << std::endl;
                    }
                }
                // std::cout << std::endl;
            }
            pcvz->printMiniZinc(toPatcher, "swarmTestFork", tFnl - tIni + 1);
            auto ppatch = pcvz->getPatch();
            if (ppatch.size() > 0) {
                pcvz->printPatch2Dot("Patch");
                p4Plan.makePatch(ppatch, 0, tIni, tFnl);
                std::cout << "patch Found" << std::endl;
                p4Plan.printTransitions();
                patchFound = true;
            } else {
                checkReach = false;
            }
        }
        // if not reachable or the cp failed, increase the patch horizon
        std::cout << "=========<Patching Iteration: Expand the Patch Horizon>==========" << std::endl;
        int ntIni = std::max(0, tIni - 1);
        int ntFnl = std::min((int) plan.size() - 1, tFnl + 1);
        if ((ntIni == tIni) && (ntFnl == tFnl)) {
            break;
        } else {
            tIni = ntIni;
            tFnl = ntFnl;
        }
    }
    // if failed, print something
    if (!patchFound) {
        std::cout << "Patch Failed" << std::endl;
    }
}

template class XSwarmTest<GR1Context, false, false>;
