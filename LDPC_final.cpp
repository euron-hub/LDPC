#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
using namespace std;

const int Z  = 96;
const int row = 12;
const int col = 24;
const int M  = row * Z; // 1152
const int N  = col * Z; // 2304
const int K  = N - M;  // 1152


int Big_H[row][col] = {
  {-1,94,73,-1,-1,-1,-1,-1,55,83,-1,-1, 7, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {-1,27,-1,-1,-1,22,79, 9,-1,-1,-1,12,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1},
  {-1,-1,-1,24,22,81,-1,33,-1,-1,-1, 0,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1},
  {61,-1,47,-1,-1,-1,-1,-1,65,25,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1},
  {-1,-1,39,-1,-1,-1,84,-1,-1,41,72,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1,-1},
  {-1,-1,-1,-1,46,40,-1,82,-1,-1,-1,79, 0,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1,-1},
  {-1,-1,95,53,-1,-1,-1,-1,-1,14,18,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1,-1},
  {-1,11,73,-1,-1,-1, 2,-1,-1,47,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1,-1},
  {12,-1,-1,-1,83,24,-1,43,-1,-1,-1,51,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1,-1},
  {-1,-1,-1,-1,-1,94,-1,59,-1,-1,70,72,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-1},
  {-1,-1, 7,65,-1,-1,-1,-1,39,49,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0},
  {43,-1,-1,-1,-1,66,-1,41,-1,-1,-1,26, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0}
};


void placeBlock(vector<vector<uint8_t>>& H, int br, int bc, int p) {
  if (p < 0) {
    return;
  } 
  const int r0 = br * Z, c0 = bc * Z;
  for (int r = 0; r < Z; ++r) {
    int c = (r + p) % Z; 
    H[r0 + r][c0 + c] = 1;
  }
}


vector<vector<uint8_t>> buildH() {
  vector<vector<uint8_t>> H(M, vector<uint8_t>(N, 0));
  for (int i = 0; i < row; ++i){
    for (int j = 0; j < col; ++j){
      placeBlock(H, i, j, Big_H[i][j]);
    }
  }
  return H;
}


bool cal_par(vector<vector<uint8_t>>& H, vector<int>& parCols) {
  parCols.assign(M, -1);
  int r = 0;
  for (int c = 0; c < N && r < M; ++c) {
    int piv = -1;
    for (int i = r; i < M; ++i) if (H[i][c]) { piv = i; break; }
    if (piv == -1) 
    continue;
    if (piv != r) swap(H[piv], H[r]);
    for (int i = 0; i < M; ++i)
      if (i != r && H[i][c])
        for (int j = c; j < N; ++j) H[i][j] ^= H[r][j];
    parCols[r] = c;
    ++r;
  }
  return (r == M);
}


vector<vector<uint8_t>> extractP(const vector<vector<uint8_t>>& H_sys,
                                 vector<int>& infoCols, const vector<int>& parCols) {
  vector<uint8_t> isParity(N, 0);
  for (int r = 0; r < M; ++r) if (parCols[r] >= 0) isParity[parCols[r]] = 1;
  infoCols.clear(); infoCols.reserve(K);
  for (int c = 0; c < N; ++c) if (!isParity[c]) infoCols.push_back(c);

  vector<vector<uint8_t>> P(M, vector<uint8_t>(K, 0));
  for (int r = 0; r < M; ++r)
    for (int j = 0; j < K; ++j)
      P[r][j] = H_sys[r][infoCols[j]];
  return P;
}


vector<uint8_t> computeParity(const vector<vector<uint8_t>>& A, const vector<uint8_t>& u) {
  vector<uint8_t> p(M, 0);
  for (int r = 0; r < M; ++r) {
    uint8_t acc = 0;
    for (int j = 0; j < K; ++j) acc ^= (A[r][j] & u[j]);
    p[r] = (acc & 1);
  }
  return p;
}


float rand49()
{  /*rand_max=7FFF (32767) */
	static int Num=0;
    double number;
    int    i;
    i=rand();
    number=(double)(i)/((unsigned) (RAND_MAX+1));
    Num++;
    if (Num >=RAND_MAX){
		time_t t;
		t=time(NULL);
		//srand((unsigned)(t%RAND_MAX));
        Num=0;
    }
    return (float)number;
}

double Normal()
{
	static int iset=0;
    static double qset;
    double vx,vy,r,temp;
    if (iset==0)//noise=normal*deviate
    {
	    do
        {
		    vx=2.0*rand49()-1.0;
            vy=2.0*rand49()-1.0;
            r =vx*vx+vy*vy;
        }while (r >=1.0 || r==0);
        temp=sqrt(-2.0*log(r)/r);
        qset=vy*temp;
        iset=1;
        return (vx*temp);
    }
    else
    {
   	    iset=0;
        return qset;
    }
}


double phi_func(double a)
{
    double result;
     double absolute;
     double tmp;
     absolute = fabs(a); 
     if( absolute < pow(10,-16)) 
     result = -37.42995; 
     else if ( absolute > 37.43) 
     result = -1.110223e-16;
     else 
     {
          tmp=(exp(a)-1)/(exp(a)+1);
          tmp=fabs(tmp);
          result= log(tmp);
     }
     return result;
}


void buildNeighbors(vector<vector<uint8_t>>& H,
                    vector<vector<int>>& CN_to_VNs,
                    vector<vector<int>>& VN_to_CNs)
{
  CN_to_VNs.assign(M, {});
  VN_to_CNs.assign(N, {});
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j)
      if (H[i][j] == 1 ) {
        CN_to_VNs[i].push_back(j);
        VN_to_CNs[j].push_back(i);
      }
    }
}


void buildPosIndex(const vector<vector<int>>& CN_to_VNs,
                   const vector<vector<int>>& VN_to_CNs,
                   vector<vector<int>>& pos_in_CNi_of_j)
{
  pos_in_CNi_of_j.assign(N, {});
  for (int j = 0; j < N; ++j) {
    pos_in_CNi_of_j[j].assign(VN_to_CNs[j].size(), -1);
  }

  for (int i = 0; i < M; ++i) {
    const auto& vlist = CN_to_VNs[i];
    for (int t = 0; t < (int)vlist.size(); ++t) {
      int j = vlist[t];
      for (int k = 0; k < (int)VN_to_CNs[j].size(); ++k) {
        if (VN_to_CNs[j][k] == i) {
          pos_in_CNi_of_j[j][k] = t;
          break;
        }
      }  
    }
  }
}


bool checkSyndrome_rows(const vector<vector<int>>& CN_to_VNs,
                        const vector<uint8_t>& bits)
{
  for (int i = 0; i < (int)CN_to_VNs.size(); ++i) {
    uint8_t acc = 0;
    const vector<int>& row = CN_to_VNs[i];   
    for (int k = 0; k < (int)row.size(); ++k) {
      int j = row[k];
      acc ^= bits[j];
    }
    if (acc & 1) {
      return false;
      }
  }
  return true;
}


struct DecodeResult {
  vector<uint8_t> hard; 
  int iters_used;
  bool parity_ok;
};


DecodeResult SPA_decode_logdomain_phi(
    const vector<vector<int>>& CN_to_VNs,
    const vector<vector<int>>& VN_to_CNs,
    const vector<vector<int>>& pos_in_CNi_of_j,
    const vector<double>& initial_LLR,
    int Iters = 20)
{
  const int Nloc = (int)VN_to_CNs.size();
  const int Mloc = (int)CN_to_VNs.size();

  vector<vector<double>> Qmsg(Nloc), Rmsg(Mloc);
  for (int j = 0; j < Nloc; ++j) {
    Qmsg[j].assign(VN_to_CNs[j].size(), initial_LLR[j]);
  }  

  for (int i = 0; i < Mloc; ++i) {
    Rmsg[i].assign(CN_to_VNs[i].size(), 0.0);
  }

  vector<double> Ltotal(Nloc, 0.0);
  vector<uint8_t> hard(Nloc, 0);

  for (int it = 0; it < Iters; ++it) {
    
    for (int i = 0; i < Mloc; ++i) {
      const auto& vlist = CN_to_VNs[i];
      int deg = (int)vlist.size();

      static thread_local vector<double> sgn, phi_pos;
      sgn.resize(deg); phi_pos.resize(deg);

      
      for (int t = 0; t < deg; ++t) {
        int j = vlist[t];
        int k = -1;
        for (int kk = 0; kk < (int)VN_to_CNs[j].size(); ++kk) {
          if (VN_to_CNs[j][kk] == i) {
            k = kk;
            break;
          }
        }  
        double q = Qmsg[j][k];
        sgn[t]     = (q >= 0.0) ? +1.0 : -1.0;
        phi_pos[t] = phi_func(q); 
      }

      
      for (int t = 0; t < deg; ++t) {
        double sign_prod = 1.0;
        double sum_phi   = 0.0;
        for (int u = 0; u < deg; ++u) {
          if (u != t) {
          sign_prod *= sgn[u];
          sum_phi   += phi_pos[u];
          }
        }  
        double mag = -phi_func(sum_phi);   
        Rmsg[i][t] = sign_prod * mag;
      }
    }

    
    for (int j = 0; j < Nloc; ++j) {
      int deg = (int)VN_to_CNs[j].size();
      double sumR = 0.0;
      for (int k = 0; k < deg; ++k) {
        int i = VN_to_CNs[j][k];
        int t = pos_in_CNi_of_j[j][k];
        sumR += Rmsg[i][t];
      }
      for (int k = 0; k < deg; ++k) {
        int i = VN_to_CNs[j][k];
        int t = pos_in_CNi_of_j[j][k];
        Qmsg[j][k] = initial_LLR[j] + (sumR - Rmsg[i][t]); 
      }
      Ltotal[j] = initial_LLR[j] + sumR;
      hard[j]   = (Ltotal[j] < 0.0) ? 1 : 0;
    }

    if (checkSyndrome_rows(CN_to_VNs, hard)) {
      return {hard, it+1, true};
    }
  }

  bool parity_ok = checkSyndrome_rows(CN_to_VNs, hard);
  return {hard, Iters, parity_ok};
}

int main() {
  srand((unsigned)time(NULL));
  double coderate = 1.0/2.0;

  auto H = buildH();
  auto H_orig = H;

  vector<int> parCols;
  if (!cal_par(H, parCols)) {
    std::cout<<"Gaussian elimination failed: rank < M.\n";
    return 1;
  }

  vector<int> infoCols;
  auto P = extractP(H, infoCols, parCols);

  
  vector<vector<int>> CN_to_VNs, VN_to_CNs, pos_in_CNi_of_j;
  buildNeighbors(H_orig, CN_to_VNs, VN_to_CNs);
  buildPosIndex(CN_to_VNs, VN_to_CNs, pos_in_CNi_of_j);

  
  mt19937_64 rng((uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count());

  std::cout<<"SNR,sets,total_bits,bit_errors,BER"<<endl;

  
  for (double SNR = 0.0; SNR <= 2.0; SNR += 0.2)
  {
    double sigma2 = 1.0 / (2.0 * coderate * pow(10.0, SNR/10.0));
    double sigma  = sqrt(sigma2);
    long long bit_errors = 0;
    int sets = 0;
    
    while (sets < 5000 || bit_errors < 100) {
      
      vector<uint8_t> u(K, 0);

      for (int i = 0; i < K; ++i) {
        u[i] = uint8_t(rng() & 1ULL);
      }  

      auto p = computeParity(P, u);
    
      vector<uint8_t> c(N, 0);

      for (int j = 0; j < K; ++j) {
        c[ infoCols[j] ] = u[j];
      }

      for (int r = 0; r < M; ++r) {
        c[ parCols[r] ]  = p[r];
      }

      vector<double> after_BPSK(N);

      for (int i = 0; i < N; ++i) {
        after_BPSK[i] = (c[i] ? +1.0 : -1.0);
      }  

      vector<double> y(N);

      for (int i = 0; i < N; ++i) {
        y[i] = after_BPSK[i] + sigma * Normal();
      }  

      vector<double> initial_LLR(N);

      for (int i = 0; i < N; ++i) {
        initial_LLR[i] = -2.0 * y[i] / sigma2;
      }  
      auto dec = SPA_decode_logdomain_phi(CN_to_VNs, VN_to_CNs, pos_in_CNi_of_j, initial_LLR, 20);
      
      for (int j = 0; j < K; ++j) {
        int bit = infoCols[j];
        if (dec.hard[bit] != u[j]) ++bit_errors;
      }
      ++sets;
    }

    long long total_bits = 1LL * sets * K;
    double BER = (total_bits==0) ? 0.0 : (double)bit_errors / (double)total_bits;

    
    std::cout<<SNR<<","<<sets<<","<<total_bits<<","<<bit_errors<<","<<BER<<endl;
  }

  return 0;
}
