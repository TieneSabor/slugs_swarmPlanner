import os
import time
import json

class regionGraph(object):
    def __init__(self):
        self.regions = ["r0","r1","r2","r3","r4"]
        self.edges = [{"r0","r1"},{"r0","r2"},{"r1","r2"},{"r1","r3"},{"r2","r4"},{"r3","r4"}]


def InitSymbols(variablesList):
    # input a variable list in json file, output two symbol dictionaries,
    # which are current and next state props
    # and show the T/F value of that symbol
    symbols = dict()
    symbolsNext = dict()
    for var in variablesList:
        symbols[var] = False
        symbolsNext[var+"'"] = False

    return symbols, symbolsNext

def safetyFunc(symbols,symbolsNext):
    # input symbols and return True or False
    s1 = not (symbols["r1"] and symbols["r4"])
    s2 = not (symbols["r0"] and symbols["r2"])

    return s1 and s2

def getPredicate(state,symbols,variablesList):
    # input state (list of 0 and 1) and symbols, output True or False
    answer = True
    for i in range(len(state)):
        if state[i] == 1:
            answer = answer and symbols[variablesList[i]]
        else:
            answer = answer and (not symbols[variablesList[i]])

    return answer


def getTrans(Nodes,ind):
    # input the node graph and the current state index,
    # output current state and next state
    curState = Nodes[str(ind)]["state"]
    nextStateInd = Nodes[str(ind)]["trans"][0]
    nextState = Nodes[str(nextStateInd)]["state"]

    return curState, nextState, nextStateInd


def getIntermediateStates(variablesList, curState, nextState, nRegion, regionGraphObj):
    # return list of pi_r and pi_r', and list of uninvolved regions

    # elements are tuples - (r,r')
    phi_r = []
    # elements are strings - r
    beta = []
    for ii in range(nRegion):
        beta.append(variablesList[ii])

    for i in range(nRegion):
        if curState[i] == 1:
            if variablesList[i] in beta:
                beta.remove(variablesList[i])
            for j in range(nRegion):
                if nextState[j] == 1:
                    if variablesList[j] in beta:
                        beta.remove(variablesList[j])
                    if {variablesList[i],variablesList[j]} in regionGraphObj.edges:
                        # NEED to check this place - conditions on False variablesList[j] at curState?
                        phi_r.append((variablesList[i],variablesList[j]))

    return phi_r, beta

def getIntermediateStatesValue(phi_r, beta, symbols, curState, nextState, variablesList):
    # return True or False given the expression of intermediate states
    phi_tau_1 = True
    phi_tau_2 = True

    for edges in phi_r:
        # print(symbols[edges[0]])
        # print(symbols[edges[1]])
        phi_tau_1 = phi_tau_1 and (symbols[edges[0]] or symbols[edges[1]])

    for regions in beta:
        phi_tau_2 = phi_tau_2 and (not symbols[regions])

    phi_q_i = getPredicate(curState,symbols,variablesList)
    phi_q_iPlus1 = getPredicate(nextState,symbols,variablesList)
    # print(phi_q_i)
    # print(phi_q_iPlus1)

    answer = phi_tau_1 and phi_tau_2 and (not phi_q_i) and (not phi_q_iPlus1)

    return answer


# def translateSafetySpec(structuredslugs_input_file_name):
    # get the safety spec from structuredslugs
    # with open(structuredslugs_input_file_name, 'r') as file:
    #     lines = file.read()
    #
    # safety = []
    # ind = lines.index("#spec: safety\n")
    # while lines[ind+1] != "\n":
    #     safety.append(lines[ind+1])
    #     ind += 1




def testSafetySpec(symbols, symbolsNext, variablesList, curState, nextState, nRegion, regionGraphObj):
    # check whether intermediate states satisfy safety specs
    phi_r, beta = getIntermediateStates(variablesList, curState, nextState, nRegion, regionGraphObj)

    N = 2**nRegion
    for i in range(N):
        for j in range(nRegion):
            symbols[variablesList[j]] = bool(int(bin(i)[2:].zfill(nRegion)[j]))

        # print(symbols)
        # print(getIntermediateStatesValue(phi_r,beta,symbols,curState,nextState,variablesList))
        # print(safetyFunc(symbols,symbolsNext))

        if getIntermediateStatesValue(phi_r,beta,symbols,curState,nextState,variablesList):
             if not safetyFunc(symbols,symbolsNext):
                 return False

    return True

def addSafety(structuredslugs_file, curState, nextState, variablesList):
    # modifies the structuredslugs file to add safety constraints
    str1 = "("
    str2 = "!("
    for i in range(len(curState)-1):
        if curState[i] == 1:
            str1 += variablesList[i] + " & "
        else:
            str1 += "!" + variablesList[i] + " & "

        if nextState[i] == 1:
            str2 += variablesList[i] + "' & "
        else:
            str2 += "!" + variablesList[i] + "' & "

    if curState[-1] == 1:
        str1 += variablesList[-1] + ")"
    else:
        str1 += "!" + variablesList[-1] + ")"

    if nextState[-1] == 1:
        str2 += variablesList[-1] + "')"
    else:
        str2 += "!" + variablesList[-1] + "')"

    str = str1 + " -> " + str2 + "\n"

    print("adding:")
    print(str)

    f = open(structuredslugs_file,"r")
    contents = f.readlines()
    f.close()

    ind = contents.index("#check intermediate states\r\n")
    contents.insert(ind+1,str+"\r\n")

    f = open(structuredslugs_file,"w")
    contents = "".join(contents)
    f.write(contents)
    f.close()





if __name__=="__main__":
    # create a region graph object
    regionGraphObj = regionGraph()
    nRegion = len(regionGraphObj.regions)

    structuredslugs_input_file_name = "demo1_coachbot_test.structuredslugs"
    slugsin_file_name = "demo1_coachbot_test.slugsin"
    automaton_file_name = "demo1_coachbot_test.txt"

    iter = 1000

    for numIter in range(iter):

        # get the slugsin file
        cmd_runStructuredSlugsParser = "tools/StructuredSlugsParser/compiler.py AAMAS2020/" + structuredslugs_input_file_name + " > AAMAS2020/" + slugsin_file_name
        os.system(cmd_runStructuredSlugsParser)
        time.sleep(0.5)

        # synthesize an automaton
        #automaton_file_name = "demo1_coachbot_test.txt"
        cmd_synthesis = "src/slugs --explicitStrategy --jsonOutput AAMAS2020/" + slugsin_file_name + " > AAMAS2020/" + automaton_file_name
        os.system(cmd_synthesis)
        time.sleep(2.0)

        # check intermediate states
        with open("AAMAS2020/"+automaton_file_name) as json_file:
            trace = json.load(json_file)

        # for var in trace["variables"]:
        #     print(var)
        #
        # curState = trace["nodes"]["0"]["state"]
        # nextStateInd = str(trace["nodes"]["0"]["trans"][0])
        # nextState = trace["nodes"][nextStateInd]["state"]
        #
        # print(curState)
        # print(nextState)

        variablesList = trace["variables"]

        symbols, symbolsNext = InitSymbols(variablesList)

        # testing
        # curState, nextState, nextStateInd = getTrans(trace["nodes"],6)
        # print(curState)
        # print(nextState)
        # phi_r,beta = getIntermediateStates(variablesList, curState, nextState, nRegion, regionGraphObj)
        # print(phi_r)
        # print(beta)
        #
        # print(testSafetySpec(symbols, symbolsNext, variablesList, curState, nextState, nRegion, regionGraphObj))

        indTmp = 0
        indSet = set()
        flag = 0
        while indTmp not in indSet:
            curState, nextState, nextStateInd = getTrans(trace["nodes"],indTmp)

            tmp = testSafetySpec(symbols, symbolsNext, variablesList, curState, nextState, nRegion, regionGraphObj)
            if tmp == False:
                addSafety("AAMAS2020/"+structuredslugs_input_file_name, curState, nextState, variablesList)
                break

            indSet.add(indTmp)
            indTmp = nextStateInd
            if nextStateInd in indSet:
                flag = 1

        if flag == 1:
            print("complete!")
            break
