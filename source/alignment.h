/* @source PolyAlign application
 * **
 * ** PolyAlign currently calculates a global alignment of two sequences
 * ** with arbitrary rational parametrs.
 * ** version 0.1.0
 * ** Please cite: Myers and Miller, CABIOS (1989)
 * **
 * ** @author Copyright (C) 2007 Brandon Barker (brandon.barker@gmail.com)
 * ** @@
 * **
 * ** This file is part of PolyAlign.
 * **
 * ** PolyAlign is free software; you can redistribute it and/or
 * ** modify it under the terms of the GNU Lesser General Public License
 * ** as published by the Free Software Foundation; either version 3 
 * ** of the License, or (at your option) any later version.
 * **
 * ** PolyAlign is distributed in the hope that it will be useful,
 * ** but WITHOUT ANY WARRANTY; without even the implied warranty of
 * ** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * ** GNU Lesser General Public License for more details.
 * **
 * ** You should have received a copy of the GNU LESSER General Public License
 * ** along with PolyAlign; if not, write to the Free Software
 * ** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * ******************************************************************************/

#include "config.h"

//Note: Find an exact or upper bound for the 
//maximum size of the edit string.
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <math.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <fstream>
#include <map>
#include <queue>

#ifdef OMP_ACTIVE
  #include <omp.h>
#endif

#include "emboss.h"

using namespace std;

namespace Para2 
{

class WorkBench;

class AlignmentSummary
{
  public:
  AlignmentSummary(WorkBench &workbench);
  ~AlignmentSummary() {delete matrix_summary;}
  unsigned long* Summary();
  string PrintSummary();
  mpq_class Score();
  map<pair<char,char>, unsigned long> *matrix_summary;
  map<char,int> *param_idx;
  unsigned long GapStart;
  unsigned long Gap;
  private:
  WorkBench *wbref;
  mpq_class score;
  bool score_set;
  
};

//General Utilities and I/O interfaces
class WorkBench
{
  public:
  string queryFileName;
  string dbFileName;
  
  AjPMatrixf FloatSubs();
  mpq_class SubMatrix(char i, char j);               //Change to N-1 for EMBOSS
  mpq_class Gap(unsigned long N) {return N>0?gapStart+(N-1)*gap:mpq_class("0");}
  WorkBench(AjPSeq query, AjPSeq dbSeq);
  ~WorkBench() {delete [] subMatrix; delete [] ASeqChar1; delete [] ASeqChar2;} 
  void ReadParams(string file_name);
  void PrintParams();
  mpq_class WriteAlign(long long* edit, AlignmentSummary &al_sum, unsigned long Len);

  AjPSeq Sequence1;
  AjPSeq Sequence2;
  char *Seq1;
  char *Seq2;
  char *ASeqChar1;
  char *ASeqChar2;
  AjPSeq ASequence1;
  AjPSeq ASequence2;
  unsigned long Seq1_len, Seq2_len;
  map<char,int> param_idx;
  mpq_class gapStart, gap;

  private:
  mpq_class *subMatrix;
};

// Note when using comparisons, the Iter variable must always be on the right side, for example:
// Bad: Iter<int> a; if(7<a) ...
// Good: Iter<int> a; if(a>=7) ...
template <class T>
class Iter
{
  public:
  T val;
  bool forward;
  Iter(){forward=true; val=0;}
  int one(){return forward?1:-1;}
  T operator=(const T& i){val=i; return val;}
  T operator=(const Iter<T>& i){val=i.val; return val;}
  T operator+(const T& o){return forward?val+o:o-val;}
  T operator-(const T& o){return forward?val-o:val+o;}
  T operator*(const T& o){return val*o;}
  T operator/(const T& o){return val/o;}
  T& operator++(){return forward?++val:--val;}
  T& operator--(){return forward?--val:++val;}
  T operator++(int){return forward?val++:val--;} // (int) = C++ magic 
  T operator--(int){return forward?val--:val++;}
  T& operator+=(const T& i){return forward?val+=i:val-=i;}
  T& operator+=(const Iter<T>& i){return forward?val+=i.val:val-=i.val;}
  T& operator-=(const T& i){return forward?val-=i:val+=i;}
  T& operator-=(const Iter<T>& i){return forward?val-=i.val:val+=i.val;}
  bool operator<=(const T& i) const {return forward?val<=i:val>=i;}
  bool operator<=(const Iter<T>& i) const {return forward?val<=i.val:val>=i.val;}
  bool operator<(const T& i) const {return forward?val<i:val>i;}
  bool operator<(const Iter<T>& i) const {return forward?val<i.val:val>i.val;}
  bool operator>=(const T& i) const {return forward?val>=i:val<=i;}
  bool operator>=(const Iter<T>& i) const {return forward?val>=i.val:val<=i.val;}
  bool operator>(const T& i) const {return forward?val>i:val<i;}
  bool operator>(const Iter<T>& i) const {return forward?val>i.val:val<i.val;}
  bool operator==(const T& i) const {return val==i;}
  bool operator==(const Iter<T>& i) const {return val==i.val;}
  bool operator!=(const T& i) const {return val!=i;}
  bool operator!=(const Iter<T>& i) const {return val!=i.val;}
  friend ostream& operator<<(ostream& os,const Iter<T> &i ){os << i.val; return os;}
};


// Used as entries in a priority queue for processing by Gotoh.
class GlobalSubProblem
{
  public:
  GlobalSubProblem(){priority=4294967295;} // This should be sufficiently large.
  GlobalSubProblem(unsigned long pri, long long As, long long Ae, long long Bs, long long Be, \
     mpq_class g1, mpq_class g2){priority=pri; Astart=As; Aend=Ae; Bstart=Bs; Bend=Be; t1=g1; t2=g2;}
  unsigned long priority;
  long long Astart, Aend, Bstart, Bend;
  mpq_class t1, t2;
  bool operator<(const GlobalSubProblem& P) const {return priority<P.priority;}
};

/* Parallel version of Myers and Miller's Global Alignment
See the paper: Optimal Alignments in Linear Space */
class GlobalAlign : public WorkBench
{
  public:
  GlobalAlign(AjPSeq query, AjPSeq dbSeq);
  ~GlobalAlign();
  void MainLoop();  //Diff is a modified version of the recursive diff in M&M; it is run from MainLoop
  void Diff(GlobalSubProblem P);
  mpq_class GotohPhase(long long Astart, long long Aend, \
    long long Bstart, long long Bend, mpq_class tn, bool forward);     

  priority_queue<GlobalSubProblem> Problems;

  long long* Edits;


  AlignmentSummary *alignment_summary;
  mpq_class score;
  
  private:
  mpq_class* CC; //CC[j] = minimum cost of conversion of Ai*to Bj
  mpq_class* DD; //DD[j] = minimum cost of conversion of Ai*to Bj /w a del
  mpq_class* RR;
  mpq_class* SS;
};


/* Standard Needleman-Wunsch method for Global Alignment
with arbitrary models and values support.
See the paper: A general method applicable to the search for similarities in the amino acid sequence of two proteins */
class NeedlemanWunsch : public WorkBench
{
  public:
  NeedlemanWunsch(AjPSeq query, AjPSeq dbSeq);
  ~NeedlemanWunsch();
  void MainLoop();  //Diff is a modified version of the recursive diff in M&M; it is run from MainLoop

  //priority_queue<BackTrace> Backtracers;

  //long long* Edits;


  AlignmentSummary* alignment_summary;
  //mpq_class score;

  private:
  int kbest;
  mpq_class* alignMatrix;
  vector<bool>  gapMatrix;
};


}
