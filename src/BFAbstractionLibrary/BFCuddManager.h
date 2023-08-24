/*
 * BFCuddManager.h
 *
 *  Created on: 06.08.2010
 *      Author: ehlers
 */

#ifndef BFCUDDMANAGER_H_
#define BFCUDDMANAGER_H_

#include <boost/utility.hpp>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <list>
#include <vector>

#include <cudd.h>

class BFBdd;
class BFBddVarCube;
class BFBddVarVector;

class BFBddManager : boost::noncopyable {
  private:
    DdManager *mgr;

  public:
    BFBddManager(unsigned int maxMemoryInMB = 3096,
                 float reorderingMaxBlowup = 1.2f);
    ~BFBddManager();

    void setAutomaticOptimisation(bool enable);
    void setReorderingMaxBlowup(float reorderingMaxBlowup);
    BFBddVarCube computeCube(const BFBdd *vars, const int *phase, int n) const;
    BFBddVarCube computeCube(const std::vector<BFBdd> &vars) const;
    BFBddVarVector computeVarVector(const std::vector<BFBdd> &vars) const;
    BFBdd readBDDFromFile(const char *filename, std::vector<BFBdd> &vars) const;
    void writeBDDToFile(const char *filename, std::string fileprefix, BFBdd bdd,
                        std::vector<BFBdd> &vars,
                        std::vector<std::string> variableNames) const;
    void writeBDDToCNFFile(const char *filename, std::string fileprefix,
                           BFBdd bdd, std::vector<BFBdd> &vars,
                           std::vector<std::string> variableNames, int *clausesInCNFIndex, int &clauseN, int &varSize) const;
    // void groupVariables(const std::vector<BFBdd> &which);
    void printStats();

    inline BFBdd constantTrue() const;
    inline BFBdd constantFalse() const;
    inline BFBdd newVariable();
    inline BFBdd multiAnd(const std::vector<BFBdd> &parts) const;
    inline BFBdd multiOr(const std::vector<BFBdd> &parts) const;

    inline DdManager *getMgr() const { return mgr; }

    friend BFBddVarVector __make_bdd_var_vect(const std::vector<BFBdd> &from,
                                              BFBddManager *mgr);
    friend class BFBdd;
    friend class BFBddVarCube;
    friend class BFBddVarVector;
};

#endif
