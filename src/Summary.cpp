/** \file Summary.cpp
 *

 * Copyright (C) 2012-2013 Ester Pantaleo


 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by


 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of


 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.


 */

#include "Summary.h"

//this function is used to sort the pair (based on likelihood) 
bool mysort(SummaryAsPair L, const SummaryAsPair M){ 
  return ((L.first)>(M.first));
}

//this function generates for each pair <key, value> in a SummaryAsMap a pair <key, <average(value), > 
//i.e., summarize the StateMap elements into pairs that can then be sorted
//based on the average value (loglikelihood) of each element of StateMap
Summary::Summary(SummaryAsMap * myMap){
  unsigned int llksize;
  double llksum;
  
  //initialize number of steps
  number_steps=0;
  
  //for each element of the map
  for (SummaryAsMap::iterator it = (*myMap).begin(); it != (*myMap).end(); ++ it){
    llksize = (it -> second).size();
    number_steps += llksize;

    //compute average llk
    llksum = 0;
    for (unsigned int i = 0; i < llksize; i++)
      llksum += (it -> second)[i];
    llksum /= llksize;

    //populate the summary in myPairs
    myPairs.push_back(std::make_pair(llksum, std::make_pair(it -> first, llksize)));
  }

  //sort the summary in myPairs
  std::sort(myPairs.begin(), myPairs.end(), mysort);
}

void Summary::print_util(std::vector<SummaryAsPair>::iterator it, ofstream & ssfile, Arguments &A){
  
  //split 101011011110000 into +,10101,10111,10000
  stringstream ss;
  ss << (it -> second).first; 
  unsigned int position = 0, i;
  if (A.pK > 0){
    ssfile << "+,"; 
    for (i=0; i<A.pK; i++){
      ssfile << ss.str().substr(position, A.lC-1) << ","; 
      position += A.lC-1;
    } 
    ssfile << " ";
  }
  
  if (A.nK > 0){
    ssfile << "-,"; 
    for (i=0; i<A.nK; i++){
      ssfile << ss.str().substr(position, A.lC-1) << ",";
      position += A.lC-1;
    }
    ssfile << " ";
  }
}

//print SummaryAsPair to file summary.txt, and print GenePredOut
void Summary::print(Arguments &A){
  double percentage;
  //print summary.txt
  ofstream ssfile;
  ssfile.open(A.summary.c_str());
  for (std::vector<SummaryAsPair>::iterator it = myPairs.begin(); it != myPairs.end(); ++ it){
    //split 101011011110000 into +,10101,10111,10000 
    print_util(it, ssfile, A);
    percentage = 100*(double)((it -> second).second)/(double)number_steps;   
    ssfile << (it -> first) << ' ' << percentage << endl; 
  }
  ssfile.close();


  //print GenePredOut
  string s ="source \"" + A.utils + "\"; ";
  //grab line 1 in file summary.txt
  s += "bitset2GenePred " + A.argv15 + "summary.txt 1 " + A.chr + " " + A.argv4 + "| sort | uniq | sort -k 3,3r  > " + A.genePredOut;
  cout << s.c_str() << endl;
  system(s.c_str());

}

