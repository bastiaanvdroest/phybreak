#include <Rcpp.h>
#include <vector>
#include <string>
using namespace Rcpp;

// Calculation of log-likelihood of sequence data
// given a phylo-transmission tree

// [[Rcpp::export(name=".likseq")]]
double likseq(IntegerVector SNPs, IntegerVector SNPfreqs,
              IntegerVector nodeparents, NumericVector nodetimes,
              double mutrate, int Nsamples) {
  int nnodes = 2*Nsamples - 1;
  int nSNPs = SNPs.size()/Nsamples;
  NumericVector likarray(nnodes * nSNPs * 4);
  NumericVector nlens(nodetimes.size());
  NumericVector SNPsums(nSNPs);
  LogicalVector routefree(nodetimes.size());
  int curnode;
  int nextnode;
  int rootnode = 0;
  double edgelen;
  double totprob;
  double result;
  
  for(int i = 0; i < Nsamples; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      likarray[(i * nSNPs + j) * 4] = 0;
      likarray[(i * nSNPs + j) * 4 + 1] = 0;
      likarray[(i * nSNPs + j) * 4 + 2] = 0;
      likarray[(i * nSNPs + j) * 4 + 3] = 0;
      if(SNPs[i*nSNPs + j] == 1) {
        likarray[(i * nSNPs + j) * 4] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 2) {
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 3) {
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 4) {
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 5) {
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 6) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 7) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 8) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 9) {
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 10) {
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 11) {
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 12) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 13) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 14) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else if(SNPs[i*nSNPs + j] == 15) {
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        
      } else {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
      }
    }
  }
  for(int i = Nsamples; i < nnodes; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      for(int k = 0; k < 4; ++k) {
        likarray[(i * nSNPs + j) * 4 + k] = 1;
      }
    }
  }
  
  for(int i = 0; i < nlens.size(); ++i) {
    if(nodeparents[i] > 0) {
      nlens[i] = nodetimes[i] - nodetimes[nodeparents[i] - 1];
    } else {
      rootnode = i;
    }
  }
  while(rootnode >= 2 * Nsamples - 1) {
    nextnode = 0;
    while(nodeparents[nextnode] - 1 != rootnode) {
      ++nextnode;
    }
    rootnode = nextnode;
  }
  for(int i = 0; i < Nsamples; ++i) {
    routefree[i] = true;
  }
  for(int i = Nsamples; i < 2 * Nsamples - 1; ++i) {
    routefree[i] = false;
  }

  for(int i = 0; i < Nsamples; ++i) {
    curnode = i;
    nextnode = nodeparents[curnode] - 1;
    edgelen = nlens[curnode];
    while(routefree[curnode] && curnode != rootnode) {
      // while(routefree[curnode] && nextnode != -1) {
        //      if(nextnode < 2*Nsamples - 1) {
        for(int j = 0; j < nSNPs; ++j) {
          totprob = likarray[(curnode * nSNPs + j) * 4];
          for(int k = 1; k < 4; ++k) {
            totprob += likarray[(curnode * nSNPs + j) * 4 + k];
          }
          for(int k = 0; k < 4; ++k) {
            likarray[(nextnode * nSNPs + j) * 4 + k] *=
              0.25 * totprob + (likarray[(curnode * nSNPs + j) * 4 + k] -
              0.25 * totprob) * exp(-mutrate * edgelen);
          }
        }
        curnode = nextnode;
        edgelen = nlens[curnode];
        nextnode = nodeparents[curnode] - 1;
//      } else {
//        edgelen += nlens[nextnode];
//        nextnode = nodeparents[nextnode] - 1;
//      }
    }
    routefree[curnode] = true;
  }
  
  for(int j = 0; j < nSNPs; ++j) {
    SNPsums[j] = 0.25*likarray[(rootnode * nSNPs + j) * 4];
    for(int k = 1; k < 4; ++k) {
      SNPsums[j] += 0.25*likarray[(rootnode * nSNPs + j) * 4 + k];
    }
    SNPsums[j] = log(SNPsums[j]) * SNPfreqs[j];
  }
  
  result = SNPsums[0];
  for(int j = 1; j < nSNPs; ++j) {
    result += SNPsums[j];
  }
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#likseq(t(curstate$d$SNP), curstate$d$SNPfr,
#         curstate$v$nodeparents, curstate$v$nodetimes,.01,200)
*/
