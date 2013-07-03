#include "Arguments.h"

//compute gamma from log_gamma                                                                        
void exp(gsl_vector *gamma, const gsl_vector *log_gamma, unsigned int size){
  for (unsigned int i = 0; i < size; i ++)
    gsl_vector_set(gamma, i, exp(gsl_vector_get(log_gamma, i)));
}

void print_my_gsl_vector(gsl_vector *v){
  for (int hh = 0;hh < v -> size; hh ++)
    cout << gsl_vector_get(v, hh) << ",";
  cout << "\n";
}

void print_my_gsl_matrix(gsl_matrix *m){
  for (int hh = 0; hh < m -> size1; hh ++)
    for (int kk = 0; kk < m -> size2; kk ++)
      cout << gsl_matrix_get(m, hh, kk) << ",";
  cout << "\n";
}

void print_to_Abuffer(gsl_vector* lambda, const double &epsilon, const double &ll, Arguments &A){
  for (unsigned int i = 0; i < A.N; i ++){
    for (unsigned int k = 0; k < A.K; k ++){
      A.buffer << setprecision(4) << gsl_vector_get(lambda, A.K * i + k) << " ";
    }
  }
  A.buffer << setprecision(4) << epsilon << " " << fixed << ll << "\n";
}

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const double &variance, unsigned int &K){
  //log likelihood normal with 0 mean 
  double result, toreturn;

  gsl_blas_ddot(log_gamma, log_gamma, &result);
  toreturn = - result;
  gsl_blas_ddot(log_lambda, log_lambda, &result);
  toreturn += result;

  return toreturn /= 2. * variance;
}

//-lambda^2 + gamma^2
double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, Arguments &A){
  //prior ratio
  double sumsq, s, tmp, tmp2, toreturn = 0.;
  for (unsigned int k = 0; k < A.K; k ++){
    sumsq = 0.;
    s = 0.;
    for (unsigned int i = 0; i < A.N; i ++){
      //*********************************
      //to speed up use stride and ddot
      //gsl_vector_const_subvector_with_stride (const gsl_vector * v, size_t offset, size_t stride, size_t n)
      tmp = gsl_vector_get(log_lambda, A.K * i + k);
      tmp2 = gsl_vector_get(log_gamma, A.K * i + k);
      sumsq += tmp * tmp - tmp2 * tmp2;
      s += tmp2 - tmp;
    }
    toreturn += (sumsq + 2. * gsl_vector_get(log_lambda_bar, k) * s) / gsl_vector_get(Var_l, k);
  }
  return 0.5 * toreturn;
}

double diff_lln(const gsl_vector* log_lambda, const gsl_vector* log_gamma, const gsl_vector* log_lambda_bar, const gsl_vector* Var_l, unsigned int &k, unsigned int &i, Arguments &A){
  //prior ratio         
  double tmp = gsl_vector_get(log_lambda, A.K * i + k);
  double tmp2 = gsl_vector_get(log_gamma, A.K * i + k);
  return  0.5 * (tmp * tmp - tmp2 * tmp2 + 2. * gsl_vector_get(log_lambda_bar, k) * (tmp2 - tmp)) / gsl_vector_get(Var_l, k);
}

double log_likelihood(const gsl_vector *lambda, const double &epsilon, const double &log_epsilon, const gsl_matrix* X, const gsl_vector* L, Arguments &A){
  double inn[A.N], toreturn = -A.LC * epsilon * A.sumC, tmp;

  for (unsigned int i = 0; i < A.N; i ++){
    tmp = 0.;
    for (unsigned int k = 0; k < A.K; k ++)
      tmp += gsl_vector_get(L, k) * gsl_vector_get(lambda, A.K * i + k);
    toreturn -= tmp * A.C[i];
    toreturn += gsl_matrix_get(X, i, 0) * (log_epsilon + A.log_C[i]);
  }
  for (unsigned int z = 1; z < A.powK; z ++){
    for (unsigned int i = 0; i < A.N; i ++)
      inn[i] = epsilon;
    for (unsigned int k = 0; k < A.K; k ++){
      if (gsl_matrix_get(A.Z, z, k) == 1){
        for (unsigned int i = 0; i < A.N; i ++)
          inn[i] += gsl_vector_get(lambda, A.K * i + k);//modify scalar product        
      }
    }
    for (unsigned int i = 0; i < A.N; i ++)
      toreturn += gsl_matrix_get(X, i, z) * (log(inn[i]) + A.log_C[i]);
  }

  return toreturn;
}

//N=1                                                                          
void init_log_lambda(gsl_vector* log_lambda, Arguments &A){
  //initialize log_lambda                                           
  rmvnorm(A.rnd, A.VAR_l, log_lambda);
}

void init_log_lambda_v(gsl_vector* log_lambda, Arguments &A){
  double V_l, sqrtc = pow(A.c, 0.5);
  gsl_vector* log_lambda_bar = gsl_vector_alloc(A.K);
  gsl_vector* Var_l = gsl_vector_alloc(A.K);

  for (unsigned int k = 0; k < A.K; k ++){
    V_l = 1. / gsl_ran_gamma(A.rnd, A.a, 1. / A.b);
    gsl_vector_set(Var_l, k, V_l);
    V_l = pow(V_l, 0.5);
    gsl_vector_set(log_lambda_bar, k, gsl_ran_gaussian(A.rnd, V_l * sqrtc));
    for (unsigned int i = 0; i < A.N; i ++)//scalar product...                  
      gsl_vector_set(log_lambda, A.K * i + k, gsl_vector_get(log_lambda_bar, k) + gsl_ran_gaussian(A.rnd, V_l));
  }

  gsl_vector_free(log_lambda_bar);
  gsl_vector_free(Var_l);
}

enum TypeOfMove{
  PChangeSS5,
  PChangeSS3,
  PAddCenter,
  PRemoveCenter,
  PSwap,
  PRecombine,
  PNull ////////////////////////////////////////////////////////////////////////////////
} Types;
 
class SS{
public:
  unsigned int c, n, n_proposed;//c = the element in the list of coordinates (pos3 or pos5); n = number of transcripts that use that splice site 
 
  SS(unsigned int &cc){
    c = cc;
    n = 0;
    n_proposed = 0;                                                         
  }
};

bool SS_less(const SS* lhs, const SS* rhs){
  return(lhs -> c < rhs -> c);
}

bool SS_equal(const SS* lhs, const SS* rhs){
  return(lhs -> c  == rhs -> c);
}

class SSpointers{
public:
  vector<SS*> v;
  unsigned int size;
};

class Center{
public:
  SSpointers ss3Pointers, ss5Pointers;
  vector<SS*> v5, v3, v5_proposed, v3_proposed;//a vector of K-dim: my_3 is the list of 3' ss used by each of the K transcripts at center c             
  unsigned int n, n_proposed;// number of transcripts that use center, total number of transcripts  
  vector<bool> is_used, is_used_proposed; // vector of size K is_used[k]==1(0) means center is used by transcript

  void clear(){
    for (unsigned int i = 0; i < ss3Pointers.size; i ++)
      ss3Pointers.v[i] -> n = 0;
    for (unsigned int i = 0; i < ss5Pointers.size; i ++)
      ss5Pointers.v[i] -> n = 0;
    v5.clear();
    v3.clear();
    n = 0;
    is_used.clear();
  }

  //define ss3Pointers.v and ss3Pointers.size 
  Center(unsigned int &i3, unsigned int &j5, vector<SS> &ss3, vector<SS> &ss5, Arguments &A){
    //define a center between A.listCoordinate[ss3[i3].c] and A.listCoordinate[ss5[j5].c]
    //define center's possible 3' and 5' ss (in the vector v), and their sizes
    vector<SS*>::iterator it;

    for (unsigned int i = 0; i <= i3; i ++){
      for (unsigned int j = j5; j < ss5.size(); j ++){
	if (A.listCoordinate[ss5[j].c] - A.listCoordinate[ss3[i3 - i].c] <= 2 * A.MAX_EX_LEN_DIV_2){
	  if (A.listCoordinate[ss5[j].c] - A.listCoordinate[ss3[i3 - i].c] > A.MIN_EX_LEN){
	    ss3Pointers.v.push_back(&(ss3[i3 - i]));
	    ss5Pointers.v.push_back(&(ss5[j]));
	  }
	}else
	  break;
      }
    }

    sort(ss3Pointers.v.begin(), ss3Pointers.v.end(), SS_less);
    it = unique(ss3Pointers.v.begin(), ss3Pointers.v.end(), SS_equal);
    ss3Pointers.v.resize( it - ss3Pointers.v.begin() );
    sort(ss5Pointers.v.begin(), ss5Pointers.v.end(), SS_less);
    it = unique(ss5Pointers.v.begin(), ss5Pointers.v.end(), SS_equal);
    ss5Pointers.v.resize(it - ss5Pointers.v.begin());
    ss3Pointers.size = (unsigned int) ss3Pointers.v.size();
    ss5Pointers.size = (unsigned int) ss5Pointers.v.size();
    n = 0;
  }

  //get the indices i3 (and i5) of the  elements in three_primes (five_primes) that have value the_3 (the_5) 
  //return 1(0) if exon is (not) compatible with center
  bool get(unsigned int &c3, unsigned int &c5, SS*& s3, SS*& s5){
    unsigned int counter = 0, p3, p5;
    for (p3 = 0; p3 < ss3Pointers.size; p3 ++){
      s3 = ss3Pointers.v[p3];
      if ((s3 -> c) == c3){
	counter = 1;
	break;
      }
    }
    for (p5 = 0; p5 < ss5Pointers.size; p5 ++){
      s5 = ss5Pointers.v[p5];
      if ((s5 -> c) == c5){
	counter += 1;
	break;
      }
    }
    if (counter == 2)
      return 1;
    return 0;
  }

  SS* draw_SS_from_posterior(SSpointers &ssPointers, Arguments &A){
    //possible_3=3,4,7,8                                                      
    //unsigned int my_proposed_ss; //chosen  SS                                       
    double q[ssPointers.size], alpha[ssPointers.size], q_cumul, u = gsl_rng_uniform(A.rnd);
    //sample q form the posterior (dirichlet)                                 
    for (unsigned int s = 0; s < ssPointers.size; s ++)
      alpha[s] = A.alphaSS + (ssPointers.v[s] -> n_proposed);
    gsl_ran_dirichlet(A.rnd, ssPointers.size, alpha, q);
    //propose SS from the posterior                                 
    q_cumul = q[0];
    for (unsigned int s = 0; s < ssPointers.size; s ++){
      if (u < q_cumul)
        return ssPointers.v[s];
      q_cumul += q[s + 1];
    }
  }

  void init_not_used(Arguments &A){
    int u = gsl_rng_uniform_int(A.rnd, ss5Pointers.size), uu = gsl_rng_uniform_int(A.rnd, ss3Pointers.size);
    is_used.push_back(0);
    v5.push_back(ss5Pointers.v[u]);
    v3.push_back(ss3Pointers.v[uu]);
  }
};

class Transcript{
public:
  unsigned int t;
  boost::dynamic_bitset<> z;
  boost::dynamic_bitset<> modifiedExon;
  //gsl_vector * number_of_anti_junctions;
  bool operator< (const Transcript &rhs) const{
    return (z < rhs.z);
  }
  
  bool operator== (const Transcript &rhs) const{
    return (z == rhs.z);
  }
  
  virtual void to_z(vector<Center> &v, unsigned int &size, Arguments &A){
    for (unsigned int c = 0; c < size; c ++){
      if (v[c].is_used[t]){//==1
        for (unsigned int cc = (v[c].v3[t] -> c); cc < (v[c].v5[t] -> c); cc ++)
          z[A.lC - cc - 2] = 1;
      }
    }
  }

  Transcript(){
  }  

  Transcript(unsigned int &my_t, Arguments &A){
    t = my_t;
    for (unsigned int i = 0; i < A.lC - 1; i++)
      z.push_back(0); 
  }
  
  Transcript(Arguments &A){
    for (unsigned int i = 0; i < A.lC - 1; i++)
      z.push_back(0);                                      
  }

  Transcript(unsigned int &my_t, vector<Center> &v, unsigned int &size, Arguments &A){
    t = my_t;
    for (unsigned int i = 0; i < A.lC - 1; i++)
      z.push_back(0);                       
    to_z(v, size, A);
  }

  int len(Arguments &A){
    int Len = 0;

    for (unsigned int c = 0; c < A.lC - 1; c ++){
      if (z[A.lC - c - 2])//==1                                  
        Len += A.L[c];
    }
    if (Len == 0) 
      return 0;
    else 
      Len -= A.RL - 1;
    //Len = max(0, Len);

    return Len;
  }    
};

class Centers{
public:
  unsigned int size, K;
  string strand;
  vector<SS> three_primes, five_primes;
  vector<Center> v;// contains the set of centers, of size "size"          
  vector<Transcript> transcripts;
  double ar_t, ar_den_t;
  vector<double> ar, ar_den; 

  void do_nothing(){}

  void clear(){
    for (unsigned int i = 0; i < size; i++)
      v[i].clear();
  }

  void print_ar(string strand, Arguments &A){
    if (K > 0){
      cout << "ACCEPTANCE RATIO FOR "<< strand << " TRANSCRIPTS = " << ar_t / ar_den_t << "\n";
      for (unsigned int i = 0; i < (unsigned int) ar.size(); i ++)
        cout << "AR " << A.my_types[i] << " = " << ar[i] / ar_den[i] << "=" << ar[i]  << "/" << ar_den[i] << "\n";
    }
  }

  void clear_ar(Arguments &A){
    if (K > 0){
      ar.assign(A.my_types.size(),0);
      ar_den.assign(A.my_types.size(),0);
      ar_t = 0.;
      ar_den_t = 0.;
    }
  } 

  void init_proposal(){
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
  void init(string my_strand, vector<unsigned int> &pos3, vector<unsigned int> &pos5, Arguments &A){
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

  bool get_c3_c5(unsigned int &genePred5SS, unsigned int &genePred3SS, unsigned int &c3, unsigned int &c5, Arguments &A){
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

  void init2(vector<unsigned int> &genePred3SS, vector<unsigned int> &genePred5SS, Arguments &A){
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
  bool check_init(string &genePredFile, Arguments &A){
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

  void init(string &genePredFile, Arguments &A){
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
};

class ProposedTranscript: public Transcript{
public:
  void to_z(vector<Center> &v, unsigned int &size, Arguments &A){
    //z must be all zeros                                                
    for (unsigned int c = 0; c < size; c ++){
      if (v[c].is_used_proposed[t])//==1               
        for (unsigned int cc = (v[c].v3_proposed[t] -> c); cc < (v[c].v5_proposed[t] -> c); cc ++)
          z[A.lC - cc - 2] = 1;
    }
  }

  ProposedTranscript(unsigned int &my_t, Centers* centers, Arguments &A):
    Transcript(my_t, A){
    to_z(centers -> v, centers -> size, A);
  }
};

class ModifiedCenter{
public:
  unsigned int t, c;

  ModifiedCenter(){}
  ModifiedCenter(unsigned int &my_c, unsigned int &my_t){
    t = my_t;
    c = my_c;
  }
};

bool is_less_than(ModifiedCenter x, ModifiedCenter y){
  return ((x.t == y.t && x.c < y.c) || x.t < y.t);
}

bool is_equal_to(ModifiedCenter x, ModifiedCenter y){
  return ((x.c == y.c) && (x.t == y.t));
}

//if element a and b are the same return 1
template<class T>
bool swap(vector<T> &v, const unsigned int &a, const unsigned int &b){
  T tmp = v[a];
  
  if (tmp == v[b])
    return 1;
  v.at(a) = v[b];
  v.at(b) = tmp;

  return 0;
}

class Move{
public:

  static vector<TypeOfMove> move_proposed;
  static vector<ModifiedCenter> modifiedCenter;
  static vector<unsigned int> t;
  static vector<SS*> s;
  static Centers* centers;
  static int counter;
  static int ar_counter; //small_world acceptance ratio
  static int small_world;

  virtual void propose(ModifiedCenter &modCenter, Arguments &A) = 0;
  
  static void clearMove(){
    move_proposed.clear();
    modifiedCenter.clear();
    t.clear();
    s.clear();
  }

  static void sortMove(){
    vector<SS*>::iterator it;
    vector<unsigned int>::iterator It;
    vector<ModifiedCenter>::iterator iT;

    sort(s.begin(), s.end(), SS_less);
    it = unique(s.begin(), s.end(), SS_equal);
    s.resize(it - s.begin());  
    sort(t.begin(), t.end());
    It = unique(t.begin(), t.end());
    t.resize(It - t.begin());
    sort(modifiedCenter.begin(), modifiedCenter.end(), is_less_than);
    iT = unique(modifiedCenter.begin(), modifiedCenter.end(), is_equal_to);
    modifiedCenter.resize(iT - modifiedCenter.begin());
  }

  void proposeChangeSS(ModifiedCenter &modCenter, SS* &my_ss, SSpointers & ss, Arguments &A){
    SS* s_current = my_ss;

    my_ss = centers -> v.at(modCenter.c).draw_SS_from_posterior(ss, A);
    /////////////////cout<<"current="<<s_current->c<<"\n";
    ////////////////cout<<"new="<<my_ss->c<<"\n";
    if (s_current != my_ss){//not accepted yet 
      (my_ss -> n_proposed) += 1;
      (s_current -> n_proposed) -= 1;
      
      s.push_back(my_ss);
      s.push_back(s_current);
      t.push_back(modCenter.t);
      modifiedCenter.push_back(modCenter);
    }
  }

  //modCenter.c, 0
  //modCenter.c, modCenter.c
  void proposeSR(unsigned int &t1, unsigned int &c1, unsigned int c2, Arguments &A){
    //if proposing a recombination of all centers you are proposing to just swap the two transcripts 
    //which is essentially a null move
    if (c1 - c2 == (Move::centers -> size))
      return;

    bool are_same, is_same;
    //sample a transcript t2 uniformly in [0,K-1] where K is the number of transcripts
    //if t2==t1 then take the K-th transcript
    unsigned int one_different = 0, t2 = gsl_rng_uniform_int(A.rnd, (centers -> K) - 1);//!!!!!!!!!

    if (t2 == t1)//!!!!!!                                                           
      t2 = (centers -> K) - 1;//!!!!!!!!!!                            
    for (unsigned int cc = c2; cc <= c1; cc ++){//if one of them is different (is_same==0) then they are different (are_same=0)
      are_same = 1;
      is_same = swap(centers -> v.at(cc).v3_proposed, t1, t2);
      are_same *= is_same;
      is_same = swap(centers -> v.at(cc).v5_proposed, t1, t2);
      are_same *= is_same; 
      is_same = swap(centers -> v.at(cc).is_used_proposed, t1, t2);
      are_same *= is_same;
      if (are_same == 0){
	one_different ++;
	modifiedCenter.push_back(ModifiedCenter(cc, t1));
	modifiedCenter.push_back(ModifiedCenter(cc, t2));
      }
    }
    if (one_different != 0){
      t.push_back(t1);
      t.push_back(t2);
    }
  }
};

vector<TypeOfMove> Move::move_proposed = vector<TypeOfMove>();
vector<ModifiedCenter> Move::modifiedCenter = vector<ModifiedCenter>();
vector<unsigned int> Move::t = vector<unsigned int>();
vector<SS*> Move::s = vector<SS*>(); 
Centers* Move::centers = NULL;

int Move::counter = 0;
int Move::ar_counter = 0;
int Move::small_world = 0;


class NullMove: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    move_proposed.push_back(PNull);
  }
};

class ProposeChangeSS5: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    proposeChangeSS(modCenter, centers -> v.at(modCenter.c).v5_proposed.at(modCenter.t), centers -> v.at(modCenter.c).ss5Pointers, A);
    centers -> ar_den[PChangeSS5] += 1;
    move_proposed.push_back(PChangeSS5);                    
  }
};

class ProposeChangeSS3: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    proposeChangeSS(modCenter, centers -> v.at(modCenter.c).v3_proposed.at(modCenter.t), centers -> v.at(modCenter.c).ss3Pointers, A);
    centers -> ar_den[PChangeSS3] += 1;
    move_proposed.push_back(PChangeSS3); 
  }
};

class ProposeAddCenter: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    SS *s5 = centers -> v.at(modCenter.c).draw_SS_from_posterior(centers -> v[modCenter.c].ss5Pointers, A);
    SS *s3 = centers -> v.at(modCenter.c).draw_SS_from_posterior(centers -> v[modCenter.c].ss3Pointers, A);
    s3 -> n_proposed += 1;
    s5 -> n_proposed += 1;
    //cout<<"s5->c "<<s5->c<<", ";
    //cout<<"s3->c "<<s3->c<<"\n";

    centers -> v[modCenter.c].n_proposed += 1;
    centers -> v[modCenter.c].is_used_proposed.at(modCenter.t) = 1;
    centers -> v.at(modCenter.c).v3_proposed.at(modCenter.t) = s3;
    centers -> v.at(modCenter.c).v5_proposed.at(modCenter.t) = s5;
    s.push_back(s3);
    s.push_back(s5);
    t.push_back(modCenter.t);
    modifiedCenter.push_back(modCenter);
    centers -> ar_den[PAddCenter] += 1;
    move_proposed.push_back(PAddCenter);
  }
};

class ProposeRemoveCenter: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    centers -> v[modCenter.c].n_proposed -= 1;
    centers -> v[modCenter.c].is_used_proposed.at(modCenter.t) = 0;
    centers -> v[modCenter.c].v3_proposed[modCenter.t] -> n_proposed -= 1;
    centers -> v[modCenter.c].v5_proposed[modCenter.t] -> n_proposed -= 1;
    s.push_back(centers -> v[modCenter.c].v3_proposed[modCenter.t]);
    s.push_back(centers -> v[modCenter.c].v5_proposed[modCenter.t]);
    t.push_back(modCenter.t);
    modifiedCenter.push_back(modCenter);
    centers -> ar_den[PRemoveCenter] += 1;
    move_proposed.push_back(PRemoveCenter);
  }
};

class ProposeRecombine: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    proposeSR(modCenter.t, modCenter.c, 0, A);
    centers -> ar_den[PRecombine] += 1;
    move_proposed.push_back(PRecombine);
  }
};

class ProposeSwap: public Move{
public:
  void propose(ModifiedCenter &modCenter, Arguments &A){
    proposeSR(modCenter.t, modCenter.c, modCenter.c, A);
    centers -> ar_den[PSwap] += 1;
    move_proposed.push_back(PSwap);                                                  
  } 
};

class Read{
public:
  unsigned int size;
  vector<unsigned int> c5, c3; //5' ends of read, 3' ends of read
  vector <double> counts;

  Read(Arguments &A){
    counts.assign(A.N, 0.);
  }

  Read(string &my_string, Arguments &A){
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
  bool is_compatible_w(Transcript &transcript, Arguments &A){
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
};

class Reads{
public:
  int size;
  vector<vector<int> > InInterval;
  vector<Read> v;
  
  void set_zero(Arguments &A){
    InInterval.assign(A.lC - 1, vector<int>());
  }

  void get(string &my_string, int &individual, Arguments &A){
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

  void print(){
    for (unsigned int i = 0; i < size; i ++){
      for (unsigned int j = 0; j < v[i].c5.size(); j ++) cout << v[i].c5[j] << ",";
      cout << "_";
      for (unsigned int j = 0; j < v[i].c3.size(); j ++) cout << v[i].c3[j] << ",";
      cout << " ";
      for (unsigned int kk = 0; kk < v[i].counts.size(); kk ++) cout<< v[i].counts[kk] << " "; 
      cout<<"\n";
    }
  }

  void print_InInterval(){
    for (unsigned int i = 0; i < InInterval.size(); i ++){
      cout << "interval " << i << ": ";
      for (unsigned int j = 0; j < InInterval[i].size(); j ++)
        cout << InInterval[i][j] << ",";
      cout << "\n";
    }
  }
};


//get the reads from the shell: pReads=positive reads (spliced reads containing a junction with a positive sense splice signal); 
//nReads=negative reads; uRead=unsigned reads (unspliced reads, that map to the reference sequence)
class allReads{
public:
  Reads pos;
  Reads neg;
  Reads uns;
  vector<int> uCounts; // number of reads that cannot map to any transcript  

  allReads (Arguments &A){   
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
    
    cout<<"pos: \n";
    pos.print();
    cout<<"neg: \n";
    neg.print();
    cout<<"uns: \n";
    uns.print();
    //scale reads with junctions
    ////////////////////////////////////////////////////////////////////////
    if (A.scale_junction_count_by != 1){
      for (unsigned int i = 0; i < A.N; i ++){
	for (unsigned int j = 0; j < pos.size; j ++)
	  pos.v[j].counts[i] *= A.scale_junction_count_by;
	for (unsigned int j = 0; j < neg.size; j ++)
	  neg.v[j].counts[i] *= A.scale_junction_count_by;
      }
    }
    ///////////////////////////////////////////////////////////////////////
    
  }
};

//check that genePrediction is either 0 or 1 or 2
class TranscriptModel{
public:
  unsigned int K;
  Centers posCenters, negCenters;
  gsl_matrix *X;
  gsl_vector *L;
  gsl_matrix *Xproposed;
  gsl_vector *Lproposed;
  vector<boost::dynamic_bitset<> > pZ, nZ, uZ, pZproposed, nZproposed, uZproposed, pZproposed2, nZproposed2, uZproposed2;
  vector<ProposedTranscript> proposedTranscripts;
   
  void print_ar(Arguments &A){
    posCenters.print_ar("+", A);
    negCenters.print_ar("-", A);
  }

  void setZ (allReads &reads, Arguments &A){
    uZ.clear();
    pZ.clear();
    nZ.clear();
    
    //(if center is positive only work with uns and pos; if center is negative only work with uns and neg
    for (unsigned int r = 0; r < reads.uns.size; r ++){
      boost::dynamic_bitset<> Z_r(A.K, 0);
      for (unsigned int k = 0; k < A.pK; k ++)
	Z_r[A.K - k - 1] = reads.uns.v[r].is_compatible_w(posCenters.transcripts[k], A);
      for (unsigned int k = A.pK; k < A.K; k ++)
	Z_r[A.K - k - 1] = reads.uns.v[r].is_compatible_w(negCenters.transcripts[k - A.pK], A);
      uZ.push_back(Z_r);
    }
    for (unsigned int r = 0; r < reads.pos.size; r ++){
      boost::dynamic_bitset<> Z_r(A.K, 0);
      for (unsigned int k = 0; k < A.pK; k ++)
	Z_r[A.K - k - 1] = reads.pos.v[r].is_compatible_w(posCenters.transcripts[k], A);
      pZ.push_back(Z_r);
    }
    for (unsigned int r = 0; r < reads.neg.size; r ++){
      boost::dynamic_bitset<> Z_r(A.K, 0);
      for (unsigned int k = A.pK; k < A.K; k ++)
	Z_r[A.K - k - 1] = reads.neg.v[r].is_compatible_w(negCenters.transcripts[k - A.pK], A);
      nZ.push_back(Z_r);
    }
  }

  void print_Z(vector<boost::dynamic_bitset<> > &aZ){
    for (unsigned int i=0;i<aZ.size();i++)
      cout<<"read "<<i<<": "<< aZ[i]<<"\n";
  } 

  void setX(gsl_matrix* my_X, Reads &reads, vector<boost::dynamic_bitset<> > &sZ, Arguments &A){
    unsigned int z;
    for (unsigned int r = 0; r < reads.size; r ++){
      z = sZ[r].to_ulong();
      for (unsigned int i = 0; i < A.N; i ++)
        gsl_matrix_set(my_X, i, z, gsl_matrix_get(my_X, i, z) + reads.v[r].counts[i]);
    }
  }

  void setX(allReads &reads, Arguments &A){
    setZ(reads,A);
    setX(X, reads.uns, uZ, A);
    setX(X, reads.pos, pZ, A);
    setX(X, reads.neg, nZ, A);
    for (unsigned int i = 0; i < A.N; i ++)
      gsl_matrix_set(X, i, 0, gsl_matrix_get(X, i, 0) + reads.uCounts[i]);
  }

  void setL(Arguments &A){
    for (unsigned int tt = 0; tt < A.pK; tt ++)
      gsl_vector_set(L, tt, posCenters.transcripts[tt].len(A));
    for (unsigned int tt = A.pK; tt < A.K; tt ++)
      gsl_vector_set(L, tt, negCenters.transcripts[tt - A.pK].len(A));
  }

  //pos: my_AK=A.K
  //neg: my_AK=A.K+A.pK
  //uns, pos, pZ, pZproposed, 
  void setZproposed(Reads &uns, Reads &my_pos, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, unsigned int &AK, Arguments &A){
    unsigned int R;
    vector<boost::dynamic_bitset<> > is_modified_u, is_modified_p;

    for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
      is_modified_u.push_back(boost::dynamic_bitset<> (uns.size, 0));
      is_modified_p.push_back(boost::dynamic_bitset<> (my_pos.size, 0));
    }

    uZproposed = uZ;
    my_pZproposed = my_pZ;

    for (unsigned int cc = 0; cc < A.lC - 1; cc ++){
      for (unsigned int r = 0; r < uns.InInterval[cc].size(); r ++){
        R = uns.InInterval[cc][r];
	for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
	  if ( !(is_modified_u[tt][R]) ){//is_modified_u[tt][R]==0 
	    if (proposedTranscripts[tt].modifiedExon[A.lC - cc - 2]){//== 1            
	      is_modified_u.at(tt)[R] = 1;                                                  
	      unsigned int k = AK + proposedTranscripts[tt].t;
	      uZproposed.at(R)[A.K - k - 1] = uns.v[R].is_compatible_w(proposedTranscripts[tt], A);
	    }
	  }
	}                                                                     
      }
      //modify Z of positive reads falling in nexon                  
      for (unsigned int r = 0; r < my_pos.InInterval[cc].size(); r ++){
        R = my_pos.InInterval[cc][r];
	for (unsigned int tt = 0; (unsigned int) tt < proposedTranscripts.size(); tt ++){
	  if ( !(is_modified_p[tt][R]) ){//is_modified_s[tt][R]==0   
	    if (proposedTranscripts[tt].modifiedExon[A.lC - cc - 2]){//== 1             
	      is_modified_p.at(tt)[R] = 1;                                       
	      unsigned int k = AK + proposedTranscripts[tt].t;
	      my_pZproposed.at(R)[A.K - k - 1] = my_pos.v[R].is_compatible_w(proposedTranscripts[tt], A);
	    }
	  }                                                           
        }
      }
    }
  }

  //choose a transcript; choose two centers c and c2; swap indicators of center c and swap indicator of center c2 (for example if center c is included remove it)       
  void setXproposed(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, Arguments &A){
    setZproposed(reads.uns, my_pos, my_pZ, my_pZproposed, AK, A);
   
    gsl_matrix_set_zero(Xproposed);
    setX(Xproposed, reads.uns, uZproposed, A);
    setX (Xproposed, my_pos, my_pZproposed, A);
    setX (Xproposed, my_neg, my_nZ, A); 
    for (unsigned int i = 0; i < A.N; i ++)
    gsl_matrix_set(Xproposed, i, 0, gsl_matrix_get(Xproposed, i, 0) + reads.uCounts[i]);
  }

  //AK = 0 or A.pK
  void setLproposed(unsigned int AK, Arguments &A){//!
    gsl_vector_memcpy(Lproposed, L);
    for (unsigned int p = 0; p < proposedTranscripts.size(); p ++){
      gsl_vector_set(Lproposed, AK + proposedTranscripts[p].t, proposedTranscripts[p].len(A));
    }
  }

  //set K, posCenters, negCenters, A.pK, A.nK, A.K, A.powK, L, X
  TranscriptModel(string &genePredFile, allReads &reads, Arguments &A){
    unsigned int a = 0;

    if (A.genePrediction==0){
      cerr << "ERROR: option genePrediction = 0 is not implemented"; exit(1);//check this!!!!
    }

    while (a == 0){
      a = 1;
      posCenters.init("+", A.pos3, A.pos5, A);        
      negCenters.init("-", A.neg3, A.neg5, A);
      a *= posCenters.check_init(genePredFile, A);
      a *= negCenters.check_init(genePredFile, A);
      if (a == 0) cout<<"Reinitializing Centers\n"; else "Centers initialized\n";
    }

    posCenters.init(genePredFile, A);
    negCenters.init(genePredFile, A);

    //cout<<"negCenters.K="<<negCenters.K<<"\n";
    //set K, A.K, A.pK, A.nK, A.powK, Z                                            
    A.set(posCenters.K, negCenters.K);
    K = A.K;
    
    //init L and X                       
    L = gsl_vector_alloc(A.K);
    X = gsl_matrix_calloc(A.N, A.powK);
    Lproposed = gsl_vector_alloc(A.K);
    Xproposed = gsl_matrix_alloc(A.N, A.powK);

    setL(A);
    setX(reads, A);

    posCenters.init_proposal();
    negCenters.init_proposal();
  }

  TranscriptModel(string &genePredFile, Centers &old_posCenters, Centers &old_negCenters, allReads &reads, Arguments &A){
    bool a = 0;

    if (A.genePrediction==0){
      cerr << "ERROR: option genePrediction = 0 is not implemented"; exit(1);
    }

    posCenters = old_posCenters;
    negCenters = old_negCenters;

    posCenters.clear();
    negCenters.clear();

    a *= posCenters.check_init(genePredFile, A);
    a *= negCenters.check_init(genePredFile, A);
    posCenters.init(genePredFile, A);
    negCenters.init(genePredFile, A);
    
    //set K, A.K, A.pK, A.nK, A.powK, Z                         
    A.set2(posCenters.K, negCenters.K);
    K = A.K;

    //init L and X                                                                      
    L = gsl_vector_alloc(A.K);
    X = gsl_matrix_calloc(A.N, A.powK);
    //Lproposed = gsl_vector_alloc(A.K);
    //Xproposed = gsl_matrix_alloc(A.N, A.powK);

    setL(A);
    setX(reads, A);
  }

  static Move *draw_move(ModifiedCenter &modifiedCenter, Arguments &A){//***********************************************
    //draw type of move//center tt    
    //only depends on the center c and on wether or not it is used by transcript t                    
    //if transcript = transcripts[t] and centers[c] = c:    
    //it originally was center.draw_type(transcript.is_using[c], A)                         
    //with probability 10% swap or recombine (if K > 1)
    Center * center = &(Move::centers -> v[modifiedCenter.c]);
    double ProbOrder2t = 0.2, p2t = 1., p_c, u;
    
    if ((Move::centers -> K) > 1)
      p2t = gsl_rng_uniform(A.rnd);
    if (p2t > ProbOrder2t){
      p_c = gsl_ran_beta(A.rnd, 1. + (center -> n_proposed), 1. + (Move::centers -> K) - (center -> n_proposed));
      u = gsl_rng_uniform(A.rnd);
      if (u > p_c){
	if ((center -> is_used_proposed[modifiedCenter.t]) == 0){
	  ///////cout<<"center -> is_used_proposed[modifiedCenter.t]) == 0\n";
	  return new NullMove;
	}else{
	  //cout<<"ProposeRemoveCenter\n";
	  return new ProposeRemoveCenter;
	}
      }else{
	if ((center -> is_used_proposed[modifiedCenter.t]) == 1){
	  if ((center -> ss3Pointers.size) + (center -> ss5Pointers.size) > 2){
	    if (gsl_rng_uniform(A.rnd) < (double) ((center -> ss3Pointers.size) - 1) / (double) ((center -> ss3Pointers.size) + (center -> ss5Pointers.size) - 2)){
	      //cout << "ProposeChangeSS3\n";
	      return new ProposeChangeSS3;
	    }else{
	      //cout << "ProposeChangeSS5\n";
	      return new ProposeChangeSS5;
	    }
	  }
	}else{
	  //cout<<"ProposeAddCenter\n";
	  return new ProposeAddCenter;
	}
      }
      //
    }*else if (p2t < 0.05)
      return new ProposeSwap;
    else
      return new ProposeRecombine;
    return new NullMove;
  }

  void reject(){
    unsigned int c, t;

    for (unsigned int i = 0; i < Move::modifiedCenter.size(); i ++){ 
      c = Move::modifiedCenter[i].c;
      t = Move::modifiedCenter[i].t;
      Move::centers -> v.at(c).n_proposed = Move::centers -> v[c].n;
      Move::centers -> v.at(c).v3_proposed.at(t) = Move::centers -> v[c].v3[t];
      Move::centers -> v.at(c).v5_proposed.at(t) = Move::centers -> v[c].v5[t];
      Move::centers -> v.at(c).is_used_proposed = Move::centers -> v[c].is_used;
    }
    for (unsigned int i = 0; i < Move::s.size(); i ++)
      Move::s[i] -> n_proposed = Move::s[i] -> n;
  }

  void accept_reject(allReads &reads, Reads &my_pos, Reads &my_neg, vector<boost::dynamic_bitset<> > &my_pZ, vector<boost::dynamic_bitset<> > &my_pZproposed, vector<boost::dynamic_bitset<> > &my_nZ, unsigned int AK, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
    unsigned int c, t;
    double alpha, llproposed;
    
    proposedTranscripts.clear();
    for (unsigned int tt = 0; tt < (Move::t).size(); tt ++){
      t = Move::t[tt];
      ProposedTranscript transcript(t, Move::centers, A);
      transcript.modifiedExon = (transcript.z ^ Move::centers -> transcripts[t].z);
      proposedTranscripts.push_back(transcript);
    }
    //set Lproposed               
    setLproposed(AK, A);
    
    for (unsigned int k = 0; k < A.K; k ++){
      if (gsl_vector_get(Lproposed, k) < 0){
	reject();
	return;
      }
    }

    setXproposed(reads, my_pos, my_neg, my_pZ, my_pZproposed, my_nZ, AK, A);    
    llproposed = log_likelihood(lambda, epsilon, log_epsilon, Xproposed, Lproposed, A);

    alpha = llproposed;
    alpha -= ll;
    alpha = min(0.,alpha);
    if (gsl_rng_uniform(A.rnd) < exp(alpha)){//accept
      if (Move::small_world == 1)
	Move::ar_counter ++;
      ll = llproposed;
      gsl_vector_swap(L, Lproposed);
      gsl_matrix_swap(X, Xproposed);
      swap(uZ, uZproposed);
      swap(my_pZ, my_pZproposed);
      for (unsigned int i = 0; i < proposedTranscripts.size(); i ++){
	t = proposedTranscripts[i].t;
	swap(Move::centers -> transcripts[t].z, proposedTranscripts[i].z);
      }

      for (unsigned int i = 0; i < Move::modifiedCenter.size(); i ++){//copy proposal
	c = Move::modifiedCenter[i].c;
	t = Move::modifiedCenter[i].t;
	Move::centers -> v.at(c).n = Move::centers -> v[c].n_proposed;
	Move::centers -> v.at(c).v3.at(t) = Move::centers -> v[c].v3_proposed[t];
	Move::centers -> v.at(c).v5.at(t) = Move::centers -> v[c].v5_proposed[t];
	Move::centers -> v.at(c).is_used = Move::centers -> v[c].is_used_proposed;
      }
      for (unsigned int i = 0; i < Move::s.size(); i ++)
	Move::s[i] -> n = Move::s[i] -> n_proposed;

      //acceptance ratios
      Move::centers -> ar_t ++;
      for (unsigned int i = 0; i < (Move::move_proposed.size()); i ++)
	Move::centers -> ar[Move::move_proposed[i]] += 1;
    }else
      reject();
  }
 
  void accept_reject_pos(allReads &reads, gsl_vector* lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
    accept_reject(reads, reads.pos, reads.neg, pZ, pZproposed, nZ, 0, lambda, epsilon, log_epsilon, ll, A);
  }

  void accept_reject_neg(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
    accept_reject(reads, reads.neg, reads.pos, nZ, nZproposed, pZ, A.pK, lambda, epsilon, log_epsilon, ll, A);
  }

  void propose_update(unsigned int &t1, allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, unsigned int &K, Arguments &A){
    Move::clearMove();
    Move::small_world = 0;
    //////////////////// draw t
    unsigned int c1;
    double u = gsl_rng_uniform(A.rnd);
    (Move::centers -> ar_den_t) += 1;
    ////////////////// draw c
    c1 = gsl_rng_uniform_int(A.rnd, Move::centers -> size);
    ////////////////// propose move
    ModifiedCenter modifiedCenter(c1, t1);
    Move *move = draw_move(modifiedCenter, A);
    move -> propose(modifiedCenter, A);
    delete move;
    /////////////////////////////////////////////////////////////////      
    if (K > 1){
      if (u > 0.7){
	unsigned int t2 = gsl_rng_uniform_int(A.rnd, K - 1);
	if (t2 == t1) 
	  t2 = K - 1;
	ModifiedCenter modifiedCenter2(c1, t2);
	Move *move2 = draw_move(modifiedCenter2, A);
	move2 -> propose(modifiedCenter2, A);
	delete move2;
      }
    }
    
    double u_sw = gsl_rng_uniform(A.rnd);
    if (u_sw > 0.7){////////////////0.9
      Move::small_world = 1;
      Move::counter ++;
      unsigned int c2;
      for (unsigned int i = 0; i < 5; i ++){
	c2 = gsl_rng_uniform_int(A.rnd, Move::centers -> size);
	ModifiedCenter modifiedCenter3(c2, t1);
	Move *move3 = draw_move(modifiedCenter3, A);
	move3 -> propose(modifiedCenter3, A);
	delete move3;
      }
    }
    Move::sortMove();
  }

  void update(allReads &reads, gsl_vector *lambda, double &epsilon, double &log_epsilon, double &ll, Arguments &A){
    unsigned int t1 = gsl_rng_uniform_int(A.rnd, A.K);

    if (t1 < A.pK ){
      Move::centers = &posCenters;
      propose_update(t1, reads, lambda, epsilon, log_epsilon, ll, A.pK, A);
      if ((Move::modifiedCenter.size()) != 0)
	accept_reject_pos(reads, lambda, epsilon, log_epsilon, ll, A);
    }else{
      t1 -= A.pK;
      Move::centers = &negCenters;
      propose_update(t1, reads, lambda, epsilon, log_epsilon, ll, A.nK, A);
      if ((Move::modifiedCenter.size()) != 0)
	accept_reject_neg(reads, lambda, epsilon, log_epsilon, ll, A);
    }
  }

  ~TranscriptModel(){
    gsl_matrix_free(Xproposed);
    gsl_vector_free(Lproposed);
    gsl_matrix_free(X);
    gsl_vector_free(L);
  }

  void print_X(){
    print_my_gsl_matrix(X); 
  }

  void print_L(){
    print_my_gsl_vector(L);
  }

  void print(string strand, Centers &centers, Arguments &A){
    if (centers.K > 0){
      //I need to make a copy of the transcripts because I'm changing their order
      ////////////////////////////////////////
      vector<Transcript> my_transcripts = centers.transcripts;
      vector<Transcript>::iterator it;
      A.buffer << strand << ",";
      sort(my_transcripts.begin(), my_transcripts.end());
      it = unique(my_transcripts.begin(), my_transcripts.end());
      my_transcripts.resize( it - my_transcripts.begin() ); 
      
      for (it = my_transcripts.begin(); it != my_transcripts.end(); ++ it)
	A.buffer << it -> z << ",";
      A.buffer << " ";
    }
  }

  void print_to_Abuffer(const double &ll, Arguments &A){
    print("+", posCenters, A);
    print("-", negCenters, A);
    A.buffer << fixed << ll << "\n";
  }
};


void update_lambda(gsl_vector* lambda, gsl_vector* log_lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A){
  double f = 1. + (double) A.N * A.c, d = A.c / f, alpha, m, s2, v, tmp, sn2, m2, llproposed;
  gsl_vector* log_gamma = gsl_vector_alloc(A.N * A.K);
  gsl_vector* gamma = gsl_vector_alloc(A.N * A.K);
  gsl_vector* log_lambda_bar = gsl_vector_alloc(A.K);
  gsl_vector* Var_l = gsl_vector_alloc(A.K);

  //draw Var_l and log_lambda_bar 
  for (unsigned int k = 0; k < A.K; k ++){
    m = 0.;
    s2 = 0;
    for (unsigned int i = 0; i < A.N; i ++)
      m += gsl_vector_get(log_lambda, A.K * i + k);
    m /= A.N;
    for (unsigned int i = 0; i < A.N; i ++)
      s2 += (gsl_vector_get(log_lambda, A.K * i + k) - m) * (gsl_vector_get(log_lambda, A.K * i + k) - m);
    m2 = m * m;
    sn2 = A.b + 0.5 * (s2 + m2 * A.N / f);
    v = 1. / gsl_ran_gamma(A.rnd, A.a + 0.5 * ((double) A.N), 1. / sn2);
    //update variance log_lambda and log_lambda_bar                                                      
    gsl_vector_set(Var_l, k, v);
    gsl_vector_set(log_lambda_bar, k, d * m * A.N + gsl_ran_gaussian(A.rnd, pow(d * v, 0.5)));
  }
  //propose log_gamma    
  //compute alpha                                                                                               
  for (unsigned int k = 0; k < A.K; k ++){
    for (unsigned int i = 0; i < A.N; i ++){
      gsl_vector_memcpy(log_gamma, log_lambda);
      gsl_vector_set(log_gamma, A.K * i + k, gsl_vector_get(log_gamma, A.K * i + k) + gsl_ran_gaussian(A.rnd, A.sd_q));
      gsl_vector_memcpy(gamma, lambda);
      gsl_vector_set(gamma, A.K * i + k, exp(gsl_vector_get(log_gamma, A.K * i + k)));

      llproposed = log_likelihood(gamma, epsilon, log_epsilon, transcriptModel.X, transcriptModel.L, A);
      alpha = llproposed - ll;
      alpha += diff_lln(log_lambda, log_gamma, log_lambda_bar, Var_l, k, i, A);
      alpha = min(0., alpha);
      if (gsl_rng_uniform(A.rnd) < exp(alpha)){
	A.ar_l += 1;//accept                                      
	gsl_vector_swap(log_lambda, log_gamma);
	gsl_vector_swap(lambda, gamma);
	ll = llproposed;
      }
    }
  }

  //Free                                                            
  gsl_vector_free(log_gamma);
  gsl_vector_free(gamma);
  gsl_vector_free(log_lambda_bar);
  gsl_vector_free(Var_l);

  return;
}

void update_epsilon (gsl_vector* lambda, double &epsilon, double &log_epsilon, TranscriptModel &transcriptModel, double &ll, Arguments &A){
  double log_delta, delta, alpha, llproposed;
  log_delta = log_epsilon + gsl_ran_gaussian(A.rnd, A.sd_q);
  delta = exp(log_delta);
  //calculate alpha                                           
  llproposed = log_likelihood(lambda, delta, log_delta, transcriptModel.X, transcriptModel.L, A);
  alpha = llproposed - ll;
  alpha += 0.5 * (log_epsilon * log_epsilon - log_delta * log_delta - 2. * A.mu_e * (log_epsilon - log_delta)) / A.var_e;
  alpha = min(0., alpha);
  if (gsl_rng_uniform(A.rnd) < exp(alpha)){
    A.ar_e ++;//accept                                  
    log_epsilon = log_delta;
    epsilon = delta;
    ll = llproposed;
  }//else reject and stay in log_epsilon        
}

int main (int argc, const char* argv[]){
  unsigned int step;
  
  Arguments A(argc, argv);
  allReads reads(A);
  double log_epsilon, epsilon, ll;
  gsl_vector* lambda;
  gsl_vector* log_lambda;

  if (A.genePrediction == 1){
    //**************************************************************************************************     
    //    INITIALIZE TRANSCRIPT MODEL, i.e., posCenters, negCenters, X and L           
    //**************************************************************************************************     
    TranscriptModel transcriptModel(A.genePredFile, reads, A);//?? which file
    //**************************************************************************************************
    //     INITIALIZE LAMBDAs AND EPSILON                          
    //**************************************************************************************************  
    //initialize lambda and epsilon
    lambda = gsl_vector_alloc(A.N * A.K);
    log_lambda = gsl_vector_alloc(A.N * A.K);
    if (A.N == 1){
      for (unsigned int i = 0; i < A.N; i ++)
	init_log_lambda(&gsl_vector_subvector(log_lambda, A.K * i, A.K).vector, A);
    }else
      init_log_lambda_v(log_lambda, A);
    exp(lambda, log_lambda, A.K * A.N);
    log_epsilon = A.mu_e + gsl_ran_gaussian(A.rnd, A.sd_e);
    epsilon = exp(log_epsilon);
  
    //print initial value///////////////////////////////////////////////////////////////////////////////
    cout<<"initial state:\n"; transcriptModel.print_to_Abuffer(ll,A);  cout << A.buffer.str();  A.buffer.str("");
    cout << "X="; transcriptModel.print_X(); cout << "L = "; transcriptModel.print_L();
    ll = log_likelihood(lambda, epsilon, log_epsilon, transcriptModel.X, transcriptModel.L, A);
    print_to_Abuffer(lambda, epsilon, ll, A); cout<<A.buffer.str(); A.buffer.str(""); 

    //*****************************************************************************************************
    //     RUN MC; print to outdata
    //*****************************************************************************************************
    string outfile = A.argv14;
    outfile += "MC_output";
    ofstream outdata;
    outdata.open(outfile.c_str());
    
    A.clear_ar();
    
    for (step = 0; step < A.MC_STEPS; step ++){
      //if (step % 1000 == 0 && step != 0) cout<<"step "<<step<<"...\n";
     
      //update transcript
      transcriptModel.update(reads, lambda, epsilon, log_epsilon, ll, A);
      update_lambda(lambda, log_lambda, epsilon, log_epsilon, transcriptModel, ll, A);
      update_epsilon(lambda, epsilon, log_epsilon, transcriptModel, ll, A);
    
      if (step % A.MC_THIN == 0){
	transcriptModel.print_to_Abuffer(ll,A);
	outdata << A.buffer.str();
	A.buffer.str("");
      }
      if (step % 10 == 0) A.outdatallk << ll << "\n";
    }
    
    //print final state ////////////////////////////////////////////////////////////////////////////////////// 
    A.print_ar_my_update(step); transcriptModel.print_ar(A);
    cout<<"ar_small_world="<<Move::ar_counter<<"/"<<Move::counter<<"="<< (double) Move::ar_counter / (double) Move::counter <<"\n";
    cout<<"last final state:\n";cout << "X="; transcriptModel.print_X(); cout << "L = "; transcriptModel.print_L();
    transcriptModel.print_to_Abuffer(ll,A); cout << A.buffer.str();A.buffer.str("");
    print_to_Abuffer(lambda, epsilon, ll, A); cout<<A.buffer.str();A.buffer.str("");
    //free                          
    outdata.close();
    
    //****************************************************************************************************
    //   SUMMARIZE_RESULTS
    //****************************************************************************************************
    cout << A.run_output2summary.c_str()<<"\n";
    system(A.run_output2summary.c_str());

    gsl_vector_free(lambda);
    gsl_vector_free(log_lambda);
  }

  //scale back to reads with junctions                                                                                                                                                  
  ////////////////////////////////////////////////////////////////////////                                                                                                      
  if (A.scale_junction_count_by != 1){
    for (unsigned int i = 0; i < A.N; i ++){
      for (unsigned int j = 0; j < reads.pos.size; j ++)
	reads.pos.v[j].counts[i] /= A.scale_junction_count_by;
      for (unsigned int j = 0; j < reads.neg.size; j ++)
	reads.neg.v[j].counts[i] /= A.scale_junction_count_by;
    }
  }
  ///////////////////////////////////////////////////////////////////////                                                                                                       
  
  TranscriptModel optimalTranscriptModel(A.genePredOut, reads, A);//define the new A.K A.pK A.nK 

  lambda = gsl_vector_alloc(A.N * A.K);
  log_lambda = gsl_vector_alloc(A.N * A.K);
  //initialize lambda and epsilon                                
  if (A.N == 1)
    init_log_lambda(log_lambda, A);
  else
    init_log_lambda_v(log_lambda, A);
  exp(lambda, log_lambda, A.K * A.N);
  log_epsilon = A.mu_e + gsl_ran_gaussian(A.rnd, A.sd_e);
  epsilon = exp(log_epsilon);
  ll = log_likelihood(lambda, epsilon, log_epsilon, optimalTranscriptModel.X, optimalTranscriptModel.L, A);

  //**************************************************************************************************     
  //    FIND OPTIMAL LAMBDAS AND EPSILON FOR FINAL STATE; print to outdata1                  
  //************************************************************************************************** 
  A.clear_ar();
 
  for (step = 0; step < A.MC_EQ + A.MC_BURNIN; step ++){
    update_lambda(lambda, log_lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, A);
    update_epsilon(lambda, epsilon, log_epsilon, optimalTranscriptModel, ll, A);
   
    //print
    if (step > A.MC_BURNIN){
      if (step % A.MC_THIN == 0){
	print_to_Abuffer(lambda, epsilon, ll, A);
	A.outdata1 << A.buffer.str();
	A.buffer.str("");
      }
    }
    if (step % 10 == 0) A.outdatallk << ll << "\n";
  }
  //print final state
  A.print_ar_my_update(step); cout<<"final state:";optimalTranscriptModel.print_to_Abuffer(ll, A);cout<<"buffer="<<A.buffer.str()<<"\n";    
  print_to_Abuffer(lambda, epsilon, ll, A);cout<<A.buffer.str();A.buffer.str("");
  cout << "X="; optimalTranscriptModel.print_X(); cout << "L = "; optimalTranscriptModel.print_L();
   
  //free
  gsl_vector_free(lambda);
  gsl_vector_free(log_lambda);
  
  return 0;
}
