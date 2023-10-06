#ifndef __EXTENSION_SWARM_TEST_HPP
#define __EXTENSION_SWARM_TEST_HPP

#include "gr1context.hpp"
#include "swarmTest2MiniZinc.hpp"
#include "swarmTestPatcher.hpp"
#include <string>

using namespace swarmTest;
/**
 * An extension that triggers that a symbolic strategy to be extracted.
 */
template <class T, bool oneStepRecovery, bool systemGoalEncoded>
class XSwarmTest : public T {
  protected:
    // New variables
    std::string outputFilename;
    convert2MiniZinc *pcvz;
    patcher p4Plan;
    BF safetySysNoRM, initState;
    std::vector<std::vector<BF>> layers2Goals;
    // specify modification
    std::string modReg;
    int modMax;
    std::string rmvEdgeFrom, rmvEdgeTo;
    // std::string reasgnFrom, reasgnTo;
    std::vector<int> reasgnRegionState;
    int reasgnStateNum, reasgnGoalNum;
    // explicit strategy
    std::vector<std::vector<int>> explicitStrategy;
    // region state assignment
    std::vector<std::vector<std::vector<int>>> rAs;

    // Inherited stuff used
    using T::addVariable;
    using T::computeVariableInformation;
    using T::determinize;
    using T::doesVariableInheritType;
    using T::initEnv;
    using T::initSys;
    using T::livenessAssumptions;
    using T::livenessGuarantees;
    using T::mgr;
    using T::postVars;
    using T::preVars;
    using T::realizable;
    using T::safetyEnv;
    using T::safetySys;
    using T::strategyDumpingData;
    using T::varCubePostOutput;
    using T::varCubePre;
    using T::variableNames;
    using T::variables;
    using T::variableTypes;
    using T::varVectorPost;
    using T::varVectorPre;
    using T::winningPositions;

    std::vector<int> counterVarNumbersPre;
    std::vector<int>
        counterVarNumbersPost;       // only used if systemGoalEncoded==True
    int goalTransitionSelectorVar;   // only used if not systemGoalEncoded==True

    XSwarmTest<T, oneStepRecovery, systemGoalEncoded>(
        std::list<std::string> &filenames)
        : T(filenames) {}

    void init(std::list<std::string> &filenames);

    BF getVarFromCVZID(bool Pre, int cvzID);

    BF regAsn2BF(std::vector<int> rA);

    // get the region map safety guarantee in BF
    BF getRMSG();

    BF parseRegionMap(std::list<std::string> &filenames);

    void updateLayers(int goalID, int layerNum, BF data);

    // update layers2Goals
    void updateL2G();

    // return the number of layer for the state "tbc" (as BF) giveb the goalID
    // return -1 if not exist
    int checkLayer(BF tbc, int goalID);

    void printStrategyBDD(BF combinedStrategy, std::string filename);

    void safetySys2IntStateRec(BF F, std::vector<std::pair<int, int>> edgePred, std::vector<std::pair<int, bool>> literal, std::vector<std::vector<std::pair<int, bool>>> &literals);

    void safetySys2IntState();

    BF safetyTransitionSetRec(BF F, std::vector<std::pair<int, int>> edgePred);

    BF safetyTransitionSet();

    void testingBF();

    void BF2CVZClauses(BF bf, whichClause wc);

    std::vector<BF> computeSymbolicStrategy();

    void computeExplicitStrategy();

    void reallocation();

    bool patchForGoal_Legacy(int goalID, std::vector<std::pair<int, int>> locP);

    bool patchForGoal(int goalID, std::vector<std::pair<int, int>> locP);

    void updateRAS();

    BF getIntermediateStateFeedback(std::vector<std::vector<int>> transition);

    std::pair<BF, int> checkReachable(BF start, BF end, BF &isFB);

    std::vector<std::vector<int>> patchForGivenHorizon(std::vector<int> iniStateAs, std::vector<int> fnlStateAs, int expLayNumb, bool final);

    std::pair<std::vector<std::vector<int>>, int> patchForGivenHorizonDoubleStrategy(std::vector<int> iniStateAs, std::vector<int> fnlStateAs, int expFirstLayNumb, int expSecondLayNumb, int goalNumb, bool final);

  public:
    void execute() {
        T::execute();
        if (realizable) {
            if (outputFilename == "") {
                throw "Internal Error.";
            } else {
                computeAndPrintSymbolicStrategy(outputFilename);
            }
        }
    }

    /**
     * @brief Compute and print out (to stdout) a symbolic strategy that is
     * winning for the system. This function requires that the realizability of
     * the specification has already been detected and that the variables
     * "strategyDumpingData" and "winningPositions" have been filled by the
     * synthesis algorithm with meaningful data.
     * @param outputStream - Where the strategy shall be printed to.
     */
    void computeAndPrintSymbolicStrategy(std::string filename);

    static GR1Context *makeInstance(std::list<std::string> &filenames) {
        return new XSwarmTest<T, oneStepRecovery, systemGoalEncoded>(filenames);
    }
};

extern template class XSwarmTest<GR1Context, false, false>;

#endif
