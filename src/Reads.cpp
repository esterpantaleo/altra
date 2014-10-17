/** \file Reads.cpp
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

#include "Reads.h"

Read::Read(Arguments &A){
  counts.assign(A.N, 0.);
}

Read::Read(string &my_string, Arguments &A){
  counts.assign(A.N, 0.);
  vector<unsigned int> my_read = i_tokenize(my_string, ",", 1);
  for (unsigned int ll = 0; ll < (unsigned int) my_read.size() - 1; ll ++){
    c5.push_back(my_read[ll]);
    c3.push_back(my_read[ll + 1]);
    ll ++;
  }
  size = (unsigned int) c5.size();
}

//transcript=001000000010110100000000011101
//read: c5=12,15   c3=13,16 0000000000001001000000000000000 at i=13 and i=14 transcript should be 0 
bool Read::is_compatible_w(Transcript &transcript, Arguments &A){
  for (unsigned int hh = 0; hh < size; hh ++){
    //check that the intron in the read corresponds to an intron in the transcript
    if (hh > 0){
      for (unsigned int i = c3[hh - 1]; i < c5[hh]; i ++){
	if (transcript.z[A.lC - i - 2] != 0)
	  return 0;
      } 
    }
    for (unsigned int i = c5[hh]; i < c3[hh]; i ++){
      if (transcript.z[A.lC - i - 2] != 1)
	return 0;
    }
  }
  
  return 1;
}
  
void Reads::set_zero(Arguments &A){
  InInterval.assign(A.lC - 1, vector<int>());
}

void Reads::get(string &my_string, int &individual, Arguments &A){
  Read read(my_string, A);
  unsigned int i = read.c5[0], r;       
  
  for (unsigned int k = 0; k < InInterval[i].size(); k ++){
    r = InInterval[i][k];
    if (read.c3 == v[r].c3 && read.c5 == v[r].c5){
      //read type is already present
      v[r].counts[individual] ++;
      return;
    }
  }
  //if read is of a new type, save the type
  size ++;
  v.push_back(read);
  v[size].counts[individual] ++;
  for (unsigned int hh = read.c5[0]; hh < read.c3[read.size - 1]; hh ++)
    InInterval[hh].push_back(size);
}

void Reads::print(){
  for (unsigned int i = 0; i < size; i ++){
    for (unsigned int j = 0; j < v[i].c5.size(); j ++) cout << v[i].c5[j] << ",";
    cout << "_";
    for (unsigned int j = 0; j < v[i].c3.size(); j ++) cout << v[i].c3[j] << ",";
    cout << " ";
    for (unsigned int kk = 0; kk < v[i].counts.size(); kk ++) cout<< v[i].counts[kk] << " "; 
    cout<<"\n";
  }
}

void Reads::print_InInterval(){
  for (unsigned int i = 0; i < InInterval.size(); i ++){
    cout << "interval " << i << ": ";
    for (unsigned int j = 0; j < InInterval[i].size(); j ++)
      cout << InInterval[i][j] << ",";
    cout << "\n";
  }
}

allReads::allReads (Arguments &A){   
  pos.InInterval.assign(A.lC - 1, vector<int>());
  neg.InInterval.assign(A.lC - 1, vector<int>());
  uns.InInterval.assign(A.lC - 1, vector<int>());
  pos.size=-1;
  neg.size=-1;
  uns.size=-1;
  
  int individual = -1;
  string my_line, my_string;
  
  while(1){ // To get you all the lines.
    cin >> my_line;
    
    if (my_line[0] == 'n'){
      cin >> my_string;
      uns.get(my_string, individual, A);
    }else{
      if (my_line[0] == '+'){
	cin >> my_string;
	pos.get(my_string, individual, A);
      }else{
	if (my_line[0] == '-'){
	  cin >> my_string;
	  neg.get(my_string, individual, A);
	}else{
	  if (my_line[0] == 'z')
	    uCounts[individual] ++;
	  else{
	    if (my_line[0] == 'm'){
	      individual++;
	      uCounts.push_back(0);
	    }
	    else{
	      if (my_line[0] == 'e')
		break;
	      else{
		cerr <<my_line<<"\n"; 
		cerr << "ERROR: invalid input string " << my_line << "\n";
	      }
	    }
	  }
	}
      }
    }
  }
  
  pos.size+=1;
  neg.size+=1;
  uns.size+=1;
  
  if (A.VERBOSE){
    cout << "printing positive, negative and unsigned read types in a format internal to altra:" << endl;
    cout<<"pos: \n";
    pos.print();
    cout<<"neg: \n";
    neg.print();
    cout<<"uns: \n";
    uns.print();
  }
  //scale reads with junctions
  if (A.scale_junction_count_by != 1){
    for (unsigned int i = 0; i < A.N; i ++){
      for (unsigned int j = 0; j < pos.size; j ++)
	pos.v[j].counts[i] *= A.scale_junction_count_by;
      for (unsigned int j = 0; j < neg.size; j ++)
	neg.v[j].counts[i] *= A.scale_junction_count_by;
    }
  }
}
