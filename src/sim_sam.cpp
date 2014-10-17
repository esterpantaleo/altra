/** \file sim_sam.cpp
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

#include "utils.h"

int main(int argc, const char* argv[]){
//Arguments:
//GenePredFile : a string with path to the genepredFile
//g_line : a string of comma separated integers specifying the lines in the genePred file
//lambda : a string of comma separated doubles specifying lambda for each element of g_line
//RL : an integer specifying read length
//M : an integer specifying the OVERHANG constraint (at least M bases are needed to map a spliced read)
//locusStart : integer
//locusEnd : integer
  
  //initialize random seed
  unsigned long int rSeed = random_seed();
  const gsl_rng_type* Tor;
  gsl_rng* rnd;
  gsl_rng_env_setup();
  Tor = gsl_rng_default;
  rnd = gsl_rng_alloc(Tor);
  gsl_rng_set(rnd, rSeed);
  
  //get arguments
  unsigned int locusStart = atoi(argv[6]), locusEnd = atoi(argv[7]), tl, RL = atoi(argv[4]), M = atoi(argv[5]), \
    N = 1, K, Y, E, line = 0, MaxEx, read, sum;
  int cc = 0, rest;
  string chr, GenePredLine, GenePredFile = argv[1], strand, s1, s2, s(RL, 'I');
  vector<string> GenePredToken;
  vector<unsigned int> g_line = i_tokenize(argv[2], ",", true), read_v; //, counter;
  vector<double> lambda = f_tokenize(argv[3], ",", true);
  ifstream ifs(GenePredFile.c_str());
  
  //check arguments
  K = (unsigned int) lambda.size();
  if (K < 1){
    cerr << "ERROR: specify a value for lambda\n";
    exit(1);
  }
  if (g_line.size() != K){
    cerr << "ERROR: wrong input, number of components of lambda should be equal to number of lines selected in the genePred file\n";
    exit(1);
  }
  
  //genePred file is in the 0-based coordinate system while sam file is in the 1-based coordinate system
  while (getline(ifs, GenePredLine)){
    line ++;
    for (unsigned int k = 0; k < K; k ++){
      tl = g_line[k];
      if (line == tl){
	cerr << "Simulating reads from transcript\n" << GenePredLine << "\n with lambda " << lambda[k] << "\n";
	
	//get transcript k from line   
	cc ++;
	GenePredToken = string_tokenize(GenePredLine, "\t ", true);
	chr = GenePredToken[1];
	strand = GenePredToken[2];
	vector<unsigned int> ExSt = i_tokenize(GenePredToken[8], ",", true);
	vector<unsigned int> ExEn = i_tokenize(GenePredToken[9], ",", true);
	E = (unsigned int) ExSt.size();
	

	//print header (only once)
	if (cc == 1)
	  cout << "@HD\tVN:1.3\tSO:sorted" << endl << "@SQ\tSN:" << chr << "\tLN:" << locusEnd+1 << endl;
        
	//define s1
	s1 = "read_name\t";
	if (strand == "+") 
	  s1 += "0\t";
	else if (strand == "-") 
	  s1 += "16\t";
	s1 += chr;
	s1 += "\t";
	s2 = "\t*\t0\t0\t";
	s2 += s;
	s2 += "\t";
	s2 += s;
	s2 += "\tNM:i:0\tNM:i:0\tNM:i:0\tNM:i:0\tNM:i:0\tNM:i:0\tNM:i:0\tNM:i:0";
        
	//for each exon j of the E exons in transcript k
	for (unsigned int j = 0; j < E; j ++){
	 
	  //generate reads with no junctions starting at position r with r in [ ExSt[j]+1,MaxEx ]
	  MaxEx = max(ExSt[j] + 1, ExEn[j] - RL + 2);
	  for (unsigned int r = ExSt[j] + 1; r < MaxEx; r ++){
	    Y = (unsigned int) gsl_ran_poisson(rnd, lambda[k]);
	    if (Y != 0){
	      for (unsigned int y = 1; y <= Y; y ++)
		cout << s1 << r << "\t16\t" << RL << "M" << s2 << "\n";
	    }
	  }
	  if (j != E - 1){
	    //generate reads with junctions
	    
            //counter.assign(E - j - 1, 0);
	    MaxEx = max(ExSt[j] + 1, ExEn[j] - RL + M + 1 + 1);
	    
	    for (unsigned int r = MaxEx; r <= ExEn[j] - M; r ++){     
	      //define read_v
	      //read_v is a vector containing the lengths of the spliced segments that form the read (example: -- -- - read_v=(2,2,1)) 
	      read_v.clear();
              rest = RL;
	      read = ExEn[j] - r + 1;
	      read_v.push_back(read);
	      rest -= read;
	      for (unsigned int jj = j + 1; jj < E; jj ++){
		if (rest <= 0) 
		  break;
		read = min(rest, (int) (ExEn[jj] - ExSt[jj]));
		read_v.push_back(read);
		rest -= read;
	      }
	      if (rest > 0) //the read would extend over the end of the transcript; the read could not be emitted by transcript k 
		break;
	      //if read overhang < M don't print read
              bool to_continue = 0;
              for (unsigned int jj = 0; jj < read_v.size(); jj ++){
		if (read_v[jj] < M){
		  to_continue = 1;
		  break;
		}
	      }
	      if (to_continue == 1)
		continue;
 

	      //print read to sam file
	      Y = (unsigned int) gsl_ran_poisson(rnd, lambda[k]);
	      if (Y != 0){
		for (unsigned int y = 1; y <= Y; y ++){
		  cout << s1 << r << "\t16\t";
		  
		  for (unsigned int jj = 0; jj < read_v.size(); jj ++){
		    cout << read_v[jj] << "M";
		    
		    if (!(read_v.size() > 1 && jj == (unsigned int) read_v.size() - 1)){
		      //counter[jj] += 1;
		      cout << ExSt[j + 1 + jj] - ExEn[j + jj] << "N";	    
		    }
		  }
		  cout << s2 << "\tXS:A:" << strand << "\n";
		}
	      }
	    }
	  }
        }
      }
    }
  }
  if (line == 0)
    cerr << "ERROR: the file " << GenePredFile << " is empty\n";

  //free
  ifs.close();
  ifs.clear();
  gsl_rng_free(rnd);
  
  return 0;
}
