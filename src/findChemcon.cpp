#include <string>
#include <iostream>
#include <list>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>

std::string pathStr;
std::string atom;
std::list<std::string> pathList;
std::list<std::pair<std::string, int>> usedBranches;
std::pair<std::string, int> newUsedBranch;
std::vector<std::string> lastPath;
std::map<std::string, std::list<std::string>> atomNeebDict;   //globAtomLabs


//
/*Equivalent to python str.split() method. Returns std::list<std::string>*/
//
std::list<std::string> splitString2List(const std::string input, const char delim) {
  std::list<std::string> inputSplit;
  std::string segment;
  std::istringstream inputISS(input);
  while(std::getline(inputISS, segment, delim)) {
    inputSplit.push_back(segment);
  }
  return inputSplit;
}

//
/*Equivalent to python str.split() method. Returns std::vector<std::string>*/
//
std::vector<std::string> splitString2Vector(const std::string input, const char delim) {
  std::vector<std::string> inputSplit;
  std::string segment;
  std::istringstream inputISS(input);
  while(std::getline(inputISS, segment, delim)) {
    inputSplit.push_back(segment);
  }
  return inputSplit;
}

void getPath() {
  //UNOPT pathStr = "";
  pathList.empty();
  std::list<std::string> passedAtoms;
  // std::cout<<"getPath"<<std::endl;
  std::string currAtom = atom;
  // std::cout<<"currAtom = " + currAtom<<std::endl;
  int steps = -1;
  // std::cout<<"steps = " + std::to_string(steps)<<std::endl;

  std::list<std::string> atomNeebs;

  while (!(find(passedAtoms.begin(), passedAtoms.end(), currAtom) != passedAtoms.end())){         //while currAtom not in passedAtoms:
    // std::cout<<"while currAtom not in passedAtoms:"<<std::endl;
    passedAtoms.push_back(currAtom);                                                                //passedAtoms.append(currAtom)
    // std::cout<<"passedAtoms.append(currAtom)"<<std::endl;                                                                                 //steps += 1
    steps++;
    // std::cout<<"steps += 1\nsteps = " + std::to_string(steps)<<std::endl;

    std::list<std::string> atomNeebs = atomNeebDict[splitString2Vector(currAtom,',')[0]];
    // std::cout<<"atomNeebs = atomNeebDict[currAtom.split(',')[0]]"<<std::endl;

    for (std::list<std::string>::const_iterator neebIter = atomNeebs.begin(), end = atomNeebs.end(); neebIter != end; ++neebIter){
      // std::cout<<"for neeb in atomNeebs:"<<std::endl;
      // std::cout<<"neeb = " + *neebIter<<std::endl;
      // UNOPT std::ostringstream branchTagOSS;
      std::pair<std::string, int> branchTag {*neebIter, steps};
      // UNOPT std::string branchTag = branchTagOSS.str();
      // std::cout<<"branchTag = " + branchTag<<std::endl;

      //if neeb not in passedAtoms and branchTag not in usedBranches:
      if ((!(find(passedAtoms.begin(), passedAtoms.end(), *neebIter) != passedAtoms.end()))
          && (!(find(usedBranches.begin(), usedBranches.end(), branchTag) != usedBranches.end()))) {
          // std::cout<<"if neeb not in passedAtoms and branchTag not in usedBranches:" + *neebIter<<std::endl;
          //
          //   std::cout<<"try:"<<std::endl;
            std::list<std::pair<std::string, int>> usedBranchesClipped;
            // std::cout<<"usedBranchesClipped"<<std::endl;
            // std::cout.flush();
            if (steps < lastPath.size()){
              if (lastPath[steps] != *neebIter){
                // std::cout<<"if lastPath.at(steps) != *neebIter"<<std::endl;
                for (std::list<std::pair<std::string, int>>::const_iterator ubIter = usedBranches.begin(), end = usedBranches.end(); ubIter != end; ++ubIter){
                  // std::cout<<"try ubIter" + *ubIter<<std::endl;
                  // std::cout<<splitString2Vector(*ubIter, '~')[0]<<std::endl;
                  // std::cout.flush();
                  if (std::get<1>(*ubIter) <= steps){
                    // std::cout<<"Clipping " + *ubIter;
                    usedBranchesClipped.push_back(*ubIter);
                  }
                }
                usedBranches = usedBranchesClipped;
              }
            }
            else {
              // std::cout<<"catch"<<std::endl;
              // std::cout.flush();

              // std::cout<<"usedBranchesClipped;";
              for (std::list<std::pair<std::string, int>>::const_iterator ubIter = usedBranches.begin(), end = usedBranches.end(); ubIter != end; ++ubIter){
                // std::cout<<"except ubIter" + *ubIter<<std::endl;
                if (std::get<1>(*ubIter) <= steps){
                  usedBranchesClipped.push_back(*ubIter);
                }
              }
            usedBranches = usedBranchesClipped;
            }

          newUsedBranch = branchTag;
          pathList.push_back(*neebIter);
          //UNOPT pathStr += *neebIter + '|';
          currAtom = *neebIter;
          break;
      }
    }
    }
  //UNOPT pathStr = pathStr.substr(0,pathStr.size()-1);
  usedBranches.push_back(newUsedBranch);
  }

//
/*Find all paths from starting atom.*/
//
std::string findAllPaths(){
  std::string paths;

  while (true) {
    getPath();
    if (pathStr.empty()) {
      break;
    }
    else {
      //UNOPT lastPath = splitString2Vector(pathStr, '|');
      lastPath.empty();
      for(std::list<std::string>::const_iterator pathIter = pathList.begin(), end = pathList.end(); pathIter != end; ++ pathIter){
        std::string tempPath = *pathIter;
        paths += tempPath.substr(0,2);
        lastPath.push_back(*pathIter);

      }
      paths += '^';
      // std::cout<<"paths" + paths<<std::endl;
    }
  }
  paths = paths.substr(0,paths.size()-1);
  return paths;
}

int main() {
  //Take raw input from python (atom and globAtomLabs dictionary)
  //Each key/value in dictionary is separated by & and neighbours are separated by *
  //i.e. "C(1)&O(1)%C(1)*C(2)" = ATOM = C(1), globAtomLabs = {O(1):[C(1),C(2)]}
  std::string input;
  std::cin>>input;
  std::list<std::string> inputSplit;                            //Raw input from python split by & into list.
  std::string segment;
  std::stringstream inputSS(input);
  int i = 0;
  std::string currAtom;
  while(std::getline(inputSS, segment, '&')) {
       inputSplit.push_back(segment);
       if (i > 0){
         if (i%2 == 0){
           std::list<std::string> neebSplit = splitString2List(segment,'*');  //Get list of neighbours for current atom.
           atomNeebDict.insert(make_pair(currAtom, neebSplit));               //Build atom neighbours map.
         }
        else {
             currAtom = segment;
           }
       }
       else {
         atom = segment;
         currAtom = atom;
       }
       i+=1;
    }

  //Testing
  // std::string neebStr;
  // for (std::list<std::string>::const_iterator iterator = atomNeebDict[atom].begin(), end = atomNeebDict[atom].end(); iterator != end; ++iterator) {
  //     neebStr += *iterator;
  //     neebStr += ", ";
  //   }
  // neebStr=neebStr.substr(0, neebStr.size()-2);
  //
  // std::cout<<"ATOM: " + atom<<std::endl;
  // std::cout<<"NEIGHBOURS: " + neebStr<<std::endl;
  //findAllPaths(atom, atomNeebDict);

  std::string x = findAllPaths();
  std::cout<<x;
  std::cout.flush();
  return 0;
}
