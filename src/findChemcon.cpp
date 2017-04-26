#include <string>
#include <iostream>
#include <list>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>

std::string pathStr;
std::string atom;
std::list<std::string> passedAtoms;
std::list<std::string> usedBranches;
std::map<std::string, std::list<std::string>> atomNeebDict;   //globAtomLabs

//
/*Equivalent to python str.split() method. Returns std::list<std::string>*/
//
std::list<std::string> splitString2List(std::string input, char delim) {
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
std::vector<std::string> splitString2Vector(std::string input, char delim) {
  std::vector<std::string> inputSplit;
  std::string segment;
  std::istringstream inputISS(input);
  while(std::getline(inputISS, segment, delim)) {
    inputSplit.push_back(segment);
  }
  return inputSplit;
}

void getPath() {

  std::string currAtom = atom;
  int steps = -1;
  std::list<std::string> atomNeebs;

  while (!(find(passedAtoms.begin(), passedAtoms.end(), currAtom) != passedAtoms.end())){         //while currAtom not in passedAtoms:

    passedAtoms.push_back(currAtom);                                                                //passedAtoms.append(currAtom)                                                                                        //steps += 1
    steps++;

    std::list<std::string> atomNeebs = atomNeebDict[splitString2Vector(currAtom,',')[0]];

    for (std::list<std::string>::const_iterator i = atomNeebs.begin(), end = atomNeebs.end(); i != end; ++i){
      std::ostringstream branchTagOSS;
      branchTagOSS<<*i<<'~'<<steps;
      std::string branchTag = branchTagOSS.str();
      std::cout << branchTag <<std::endl;
    }
    }

  }



std::vector<std::string> lastPath;


//
/*Find all paths from starting atom.*/
//
std::string findAllPaths(std::string atom, std::map<std::string, std::list<std::string>> atomNeebDict) {


  while (true) {

  }
  getPath();
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
  getPath();
  std::cout.flush();
  return 0;
}
