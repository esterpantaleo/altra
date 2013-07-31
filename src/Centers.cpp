#include "Centers.h"

void Centers::clear(){
  for (unsigned int i = 0; i < size; i++)
    v[i].clear();
}

void Centers::print_ar(string strand, Arguments &A){
  if (K > 0){
    cout << "ACCEPTANCE RATIO FOR "<< strand << " TRANSCRIPTS = " << ar_t / ar_den_t << "\n";
    for (unsigned int i = 0; i < (unsigned int) ar.size(); i ++)
      cout << "AR " << A.my_types[i] << " = " << ar[i] / ar_den[i] << "=" << ar[i]  << "/" << ar_den[i] << "\n";
  }
}

void Centers::clear_ar(Arguments &A){
  if (K > 0){
    ar.assign(A.my_types.size(),0);
    ar_den.assign(A.my_types.size(),0);
    ar_t = 0.;
    ar_den_t = 0.;
  }
} 

void Centers::init_proposal(){
  if (K > 0){
      for (unsigned int s = 0; s < three_primes.size(); s ++)
	three_primes[s].n_proposed = three_primes[s].n;
      for (unsigned int s = 0; s < five_primes.size(); s ++)
	five_primes[s].n_proposed = five_primes[s].n;
      for (unsigned int c = 0; c < v.size(); c++){	
	v[c].v5_proposed = v[c].v5;
	v[c].v3_proposed = v[c].v3;
	v[c].n_proposed = v[c].n;
	v[c].is_used_proposed = v[c].is_used;
      }
  }
}

//set strand, three_primes, five_primes, v, size  
void Centers::init(string my_strand, vector<unsigned int> &pos3, vector<unsigned int> &pos5, Arguments &A){
  unsigned int p3, p5;
  strand = my_strand;
  three_primes.clear();
  five_primes.clear();
  v.clear();
  size = 0;
  
  if (!(pos3.empty() || pos5.empty())){
    //example: pos3 = 4,5,6,9 pos5=2,7,8,10
    //listCoordinate=1100000.5,1100025.5,1100050.5,110051.5                  
    //1100001,1100051 -> pos3 = 0,2 (,3,4,8,9)                                   
    //1100025,1100051 -> pos5 = 1,3 (,10,11)
    for (p3 = 0; p3 < pos3.size(); p3 ++)
      three_primes.push_back(SS(pos3[p3]));
    for (p5 = 0; p5 < pos5.size(); p5 ++)
      five_primes.push_back(SS(pos5[p5]));
    for (unsigned int c = 0; c < A.lC - 1; c ++){
      for (p3 = 0; p3 < pos3.size() - 1; p3 ++){
	if (pos3[p3] == c){
	  for (p5 = 0; p5 < pos5.size(); p5 ++){
	    if (pos5[p5] >= c && pos5[p5] < pos3[p3 + 1]){
	      Center center(p3, p5, three_primes, five_primes, A);
	      if (!(center.ss3Pointers.size == 0 || center.ss5Pointers.size == 0)){
		v.push_back(center);
		goto stop;
	      }
	    }
	  }
	}
      }
      if (pos3[p3] == c){
	for (p5 = 0; p5 < pos5.size(); p5 ++){
	  if (pos5[p5] >= c){
	    Center center(p3, p5, three_primes, five_primes, A);// set center c, i.e.:    v[c].ss3Pointers.v, v[c].ss5Pointers.v, v[c].ss3Pointers.size, v[c].ss5Pointers.size, v[c].n
	    if (!(center.ss3Pointers.size == 0 || center.ss5Pointers.size == 0)){
	      v.push_back(center);
	      goto stop;
	    }
	  }
	}
      }
    stop:
      do_nothing();
    }
  }
  size = (unsigned int) v.size();
}

bool Centers::get_c3_c5(unsigned int &genePred5SS, unsigned int &genePred3SS, unsigned int &c3, unsigned int &c5, Arguments &A){
  while (genePred5SS - genePred3SS - 1 > 2 * A.MAX_EX_LEN_DIV_2){
    cerr << "WARNING: MAX EXON LENGTH=" << A.MAX_EX_LEN_DIV_2 * 2 << " is too small to accommodate gene prediction in the genepred file.\n";
    A.MAX_EX_LEN_DIV_2 =  (genePred5SS - genePred3SS - 1)/2 + 200;
    cerr << "Setting MAX EXON LENGTH to " << 2 * A.MAX_EX_LEN_DIV_2 <<".\n"; 
    return 1;
  }
  //check that exon boundaries in genepred file are in listCoordinate          
  for (c3 = 0; c3 < A.lC; c3 ++){
    if (genePred3SS == A.listCoordinate[c3] - 0.5){
      break;
    }
  }
  if (c3 == A.lC){
    for (c3 = 0; c3 < A.lC; c3 ++){
      if (abs(genePred3SS - A.listCoordinate[c3] + 0.5) <= 3){
	break;
      }
    }
  }
  if (c3 == A.lC){
    cerr << "ERROR: exon boundary " << genePred3SS << " in genepred file is not in listCoordinate\n"; 
    exit(1);
  }
  for (c5 = c3 + 1; c5 < A.lC; c5 ++){
    if (genePred5SS == A.listCoordinate[c5] - 0.5)
      break;
  }
  if (c5 == A.lC){
    for (c5 = c3 + 1; c5 < A.lC; c5 ++){
      if (abs(genePred5SS - A.listCoordinate[c5] + 0.5) <= 3)
	break;
    }
  }
  if (c5 == A.lC){
    cerr << "ERROR: exon boundary " << genePred5SS << " in genepred file is not in listCoordinate\n"; 
    exit(1);
  }
  return 0;
}

void Centers::init2(vector<unsigned int> &genePred3SS, vector<unsigned int> &genePred5SS, Arguments &A){
  unsigned int c3, c5, c, next_c = 0;
  SS* s3;
  SS* s5;
  
  for (unsigned int e = 0; e < (unsigned int) genePred3SS.size(); e ++){//for each exon
    get_c3_c5(genePred5SS[e], genePred3SS[e], c3, c5, A);
    for (c = next_c; c < size; c++){
      if (v[c].get(c3, c5, s3, s5)){
	next_c = c + 1;
	v[c].is_used.push_back(1);
	v[c].n += 1;
	v[c].v5.push_back(s5);
	v[c].v3.push_back(s3);
	(s3 -> n) += 1;
	(s5 -> n) += 1;
	break;
      }else//init from prior
	v.at(c).init_not_used(A);
    }
    if (c == size){
      cerr <<"ERROR: exon "<< e <<" in genepred file does not correspond to any center\n 3: "<<genePred3SS[e]<<", 5: " << genePred5SS[e]<<"\n";
      exit(1);
    }
  }
  for (unsigned int cc = max(c, next_c); cc < size; cc ++)
    v.at(cc).init_not_used(A);
}

//set K, for c in 1,...,size set v[c].v5, v[c].v3, v[c].is_used and the n-s
bool Centers::check_init(string &genePredFile, Arguments &A){
  //set K (A.pK = number of pos transcripts in genePredFile, A.nK = number of negative transcripts in genePredFile )     
  bool a = 1;
  string genePredLine;
  vector<unsigned int> genePred3SS, genePred5SS;
  vector<string> genePredToken;
  
  K = 0;
  if (size == 0)
    return 1;
  ifstream ifs(genePredFile.c_str());
  while (getline(ifs, genePredLine)){
    genePredToken = string_tokenize(genePredLine, "\t ", true);
    if (genePredToken[2].c_str() == strand){
      genePred3SS = i_tokenize(genePredToken[8], ",", 1);//3' SS of the gene in File                    
      genePred5SS = i_tokenize(genePredToken[9], ",", 1);//5' SS of the gene in File 
      if (genePred3SS.size() != genePred5SS.size()){
	cerr <<"ERROR: genePred file contains an incomplite list of splice sites\n";
	exit(1);
      }
      for (unsigned int hh = 0; hh < genePred3SS.size(); hh ++){
	while (genePred5SS[hh] - genePred3SS[hh] - 1 > 2 * A.MAX_EX_LEN_DIV_2){
	  cerr << "WARNING: MAX EXON LENGTH=" << A.MAX_EX_LEN_DIV_2 * 2 << " is too small to accommodate gene prediction in the genepred file.\n";
	  A.MAX_EX_LEN_DIV_2 =  (genePred5SS[hh] - genePred3SS[hh] - 1)/2 + 200;
	  cerr << "Setting MAX EXON LENGTH to " << 2 * A.MAX_EX_LEN_DIV_2 <<".\n";
	  a *= 0;
	}
      }
      K += 1;
    }
  }
  if (a == 0)
    return 0;
  ifs.close();
  ifs.clear();
  
  return 1;
}

void Centers::init(string &genePredFile, Arguments &A){
  string genePredLine;
  vector<unsigned int> genePred3SS, genePred5SS;
  vector<string> genePredToken;
  unsigned int t;
  
  ifstream ifs2(genePredFile.c_str());
  if (K != 0){
    t = 0;
    while (getline(ifs2, genePredLine)){
      //cout<<genePredLine<<"\n";
      genePredToken = string_tokenize(genePredLine, "\t ", true);
      if (genePredToken[2].c_str() == strand){
	genePred3SS = i_tokenize(genePredToken[8], ",", 1);//3' SS of the gene in File  
	genePred5SS = i_tokenize(genePredToken[9], ",", 1);//5' SS of the gene in File
	//initialize the centers at the t-h transcript                     
	init2(genePred3SS, genePred5SS, A);
	transcripts.push_back(Transcript(t, v, size, A));//??
	t ++;
      }
    }
  }
  ifs2.close();
  ifs2.clear();
  clear_ar(A);
}
