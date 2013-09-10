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


#include "alignment.h" 
using namespace std;
using namespace Para2;

AlignmentSummary::AlignmentSummary(WorkBench &workbench)
{
  GapStart=0;
  Gap=0;
  wbref = &workbench; 
  param_idx=&(wbref->param_idx);
  matrix_summary = new map<pair<char,char>, unsigned long>;
  map<char,int>::iterator i;
  map<char,int>::iterator j;
  for(i=param_idx->begin();i!=param_idx->end();i++)
    for(j = param_idx->begin();j!=param_idx->end();j++)
      matrix_summary->insert(make_pair(make_pair(i->first,j->first),0));
}

/* Order of output will always be:
  (GapStart,Gaps,H_i,j) as i->N and j->N.
  where i is row# and j is col# and N is matrix length. */

unsigned long* AlignmentSummary::Summary()
{
  score_set = false;
  unsigned long *summary;
  summary = new unsigned long[2+matrix_summary->size()];
  summary[0] = GapStart;
  summary[1] = Gap;
  map<char,int>::iterator i;
  map<char,int>::iterator j;
  unsigned long k=2;
  for(i = param_idx->begin();i!=param_idx->end();i++)
    for(j = param_idx->begin();j!=param_idx->end();j++)
    {
      summary[k]=(*matrix_summary)[make_pair(i->first,j->first)];  
      k++; 
    }
  return summary;  
}


mpq_class AlignmentSummary::Score()
{
  if(!score_set)
  {
    score=0;
    map<char,int>::iterator i;
    map<char,int>::iterator j;
    for(i = param_idx->begin();i!=param_idx->end();i++)
      for(j = param_idx->begin();j!=param_idx->end();j++)
      {
        score += wbref->SubMatrix(i->first,j->first) * (*matrix_summary)[make_pair(i->first,j->first)];
      }
    score += wbref->gapStart * GapStart;
    score += wbref->gap * Gap;

    score_set=true;
  }
  
  return score;

}


string AlignmentSummary::PrintSummary()
{
  unsigned long *summary;
  string *param_names;
  unsigned long params_size;
  string summary_text="";

  params_size=2+matrix_summary->size();
  cout << "matrix_summary->size: " << matrix_summary->size() << endl;
  cout << "params_size: " << params_size << endl;
  summary = new unsigned long[params_size];
  param_names = new string[params_size];
  summary[0] = GapStart;
  summary[1] = Gap;
  param_names[0] = "GapStart";
  param_names[1] = "Gap";
  map<char,int>::iterator i;
  map<char,int>::iterator j;
  unsigned long k=2;
  for(i = param_idx->begin();i!=param_idx->end();i++)
    for(j = param_idx->begin();j!=param_idx->end();j++)
    {
      char tmp[2]; //why 2???????????????????????????????
      summary[k]=(*matrix_summary)[make_pair(i->first,j->first)];
      tmp[0]=i->first;
      tmp[1]=j->first;
      tmp[2]='\0';
      param_names[k] = tmp;
      k++;
    }
  for(unsigned long i=0;i<params_size;i++)
  {
    char tmp[25]; //Fix this to be general as opposed to arbitrary.
    sprintf(tmp,"%lu",summary[i]);
    summary_text += (param_names[i] + ":\t" + tmp + "\n");
  }
  
  delete [] param_names;
  delete [] summary;
  return summary_text;
}



GlobalAlign::GlobalAlign(AjPSeq query, AjPSeq dbSeq):WorkBench(query, dbSeq)
{
  ReadParams("params.txt"); //Change this to a command line option later.
  CC = new mpq_class[Seq2_len+2];
  DD = new mpq_class[Seq2_len+2];
  RR = new mpq_class[Seq2_len+2];
  SS = new mpq_class[Seq2_len+2];
  Edits = new long long[Seq1_len+Seq2_len];
  for (unsigned long i=0;i<Seq1_len+Seq2_len;i++)
    Edits[i]=-1;

  CC[0]="0"; 
  RR[Seq2_len+1]="0";
  alignment_summary = new AlignmentSummary(*this);
}

GlobalAlign::~GlobalAlign()
{
  delete [] CC;
  delete [] DD;
  delete [] RR;
  delete [] SS;
  delete alignment_summary;
  delete [] Edits;
}

mpq_class GlobalAlign::GotohPhase(long long Astart, long long Aend, \
  long long Bstart, long long Bend, mpq_class tn, bool forward)
{
  mpq_class* C;
  mpq_class* D;
  mpq_class C0;
  mpq_class sS_ij;
  C0 = "0";
  Iter<long long> i, j;
  if(forward)
  {
    C=CC;
    D=DD;
  }
  else
  {
    C=RR;
    D=SS;
    i.forward=false;
    j.forward=false;
  }
  mpq_class e,c,s,t;
  t = -(gapStart-gap); //Don't charge first gap count
  for(j=Bstart;j<=Bend;j++)
  {
    C[j.val]=t=t-gap;
    D[j.val]=t-gapStart;
  } 
  t=tn;
  for(i=Astart;i<=Aend;i++)
  {
    
    s=C0;
    C0=c=t=t-gap;
    e=t-gapStart;
    for(j=Bstart;j<=Bend;j++)              
    {                                      
      sS_ij= s-SubMatrix(Seq1[i.val-1], Seq2[j.val-1]);
      e = ((e-gap) < c-gapStart ? (e-gap):(c-gapStart)); //Gap(N)=h*(N-1)+g
      D[j.val]=((D[j.val]-gap) < (C[j.val]-gapStart) ? (D[j.val]-gap):(C[j.val]-gapStart)); //Gap(N)=h*(N-1)+g  
      if (D[j.val] <= e && D[j.val] <= sS_ij)  
        c=D[j.val];
      else if (e <= D[j.val] && e <= sS_ij)
        c=e;
      else
      {
        c=sS_ij;
      }
      s=C[j.val];
      C[j.val]=c;
    }
  }
  return C0;
}

void GlobalAlign::MainLoop()
{
  //Create a Prioriy Queue to store sub problems in.
  vector<GlobalSubProblem> ActiveQueue;
  GlobalSubProblem initialProblem(0,1,Seq1_len,1,Seq2_len,-(gapStart-gap),-(gapStart-gap));
  Problems.push(initialProblem);
  #ifdef OMP_ACTIVE
  int nt = omp_get_max_threads();
  #endif
  while(!Problems.empty())
  {
    unsigned long P_size = Problems.size();
    //Separate queue used to reduce critical section time.
    ActiveQueue.resize(P_size);
    for(unsigned long i=0;i<P_size;i++)
    {
      ActiveQueue[i]=Problems.top();
      Problems.pop();
    }
    
    #pragma omp parallel for num_threads(nt/2)
    for(unsigned long i=0;i<P_size;i++)
    {
      Diff(ActiveQueue[i]);
    }
  }

}

void GlobalAlign::Diff(GlobalSubProblem P)  
{
  long long split = (P.Aend+P.Astart-1)/2;
  Iter<long long> i, j;
  mpq_class CC0, RR0;
  j.forward = true; i.forward=true; //Necessary?
  if (P.Bend < P.Bstart) 
  {
    if (P.Aend >= P.Astart)  
    {
      Edits[split+P.Bend]=P.Aend-P.Astart+2; 
    }
  }
  else if (P.Aend < P.Astart) 
  {
    Edits[split+P.Bend]=-(P.Bend-P.Bstart+1)-1;
  }
  else if (P.Aend-P.Astart == 0) 
  {
    mpq_class min_B_delA, min_insB_H_insB, temp; 
    long long j_hit=P.Bstart;
    min_B_delA = (P.t1<P.t2?P.t1:P.t2)-gap-Gap(P.Bend-P.Bstart+1); 
    min_insB_H_insB =  /* Gap(0) */ - SubMatrix(Seq1[P.Astart-1],Seq2[P.Bstart-1]) - Gap(P.Bend-P.Bstart);
    for(j=P.Bstart+1;j<=P.Bend;j++)
    {
      temp = - Gap(j.val-P.Bstart) - SubMatrix(Seq1[P.Astart-1],Seq2[j.val-1]) - Gap(P.Bend-j.val);
      if (temp <= min_insB_H_insB)
      {
        min_insB_H_insB=temp;
        j_hit = j.val;
      }
    }
    if (min_insB_H_insB<min_B_delA) 
    {
      Edits[split+j_hit-2]=-(j_hit-P.Bstart)-1;
      Edits[split+j_hit]=-(P.Bend-j_hit)-1;
      Edits[split+j_hit-1]=1;      
    }
    else
    {
      Edits[split+P.Bstart-1]=2;
      Edits[split+P.Bstart]=-(P.Bend-P.Bstart+1)-1;
    }
  }
  else
  {
    #pragma omp parallel num_threads(2)
    {
    #pragma omp sections
    {
      #pragma omp section
      {
        CC0 = GotohPhase(P.Astart,split,P.Bstart,P.Bend,P.t1,true);
      }
      #pragma omp section
      {
        RR0 = GotohPhase(P.Aend,split+1,P.Bend,P.Bstart,P.t2,false);
      }
    }
    #pragma omp barrier
    }
    long long minj;
    bool type1;
    mpq_class minT1, minT2, minCost;
    minT1=CC0+RR[P.Bstart];
    minT2=CC0+SS[P.Bstart]+gapStart-gap;
    minj = P.Bstart-1;
    if (minT1 <= minT2)
    {
      type1 = true;
      minCost = minT1;
    }
    else
    {
      type1 = false;
      minCost = minT2;
    }    //Now find min i,j
    for(j=P.Bstart;j<P.Bend;j++)
    {
      minT1=CC[j.val]+RR[j.val+1];
      minT2=DD[j.val]+SS[j.val+1]+gapStart-gap;
      if (minT1 <= minCost)
      {
        type1 = true;
        minCost = minT1;
        minj = j.val;
      }
      if (minT2 <= minCost) 
      {
        type1 = false;
        minCost = minT2;
        minj = j.val;
      }
    }
    minT1=CC[P.Bend]+RR0;
    minT2=DD[P.Bend]+RR0+gapStart-gap;
    if (minT1 <= minCost)
    {
      type1 = true;
      minCost = minT1;
      minj = j.val;
    }
    if (minT2 <= minCost)
    {
      type1 = false;
      minCost = minT2;
      minj = j.val;
    }

    if (type1)
    {
      GlobalSubProblem P1(P.priority+1,P.Astart,split,P.Bstart,minj,P.t1,-(gapStart-gap)); 
      GlobalSubProblem P2(P.priority+1,split+1,P.Aend,minj+1,P.Bend,-(gapStart-gap),P.t2);
      #pragma omp critical(enqueue_problems)
      {
        Problems.push(P1);
        Problems.push(P2);   
      }
    }
    else
    {
      GlobalSubProblem P1(P.priority+1,P.Astart,split-1,P.Bstart,minj,P.t1,0);
      GlobalSubProblem P2(P.priority+1,split+2,P.Aend,minj+1,P.Bend,0,P.t2);

      #pragma omp critical(enqueue_problems)
      {
        Problems.push(P1);
        Problems.push(P2);
      }
      Edits[split+minj]=3; 
    }
  }
}


mpq_class WorkBench::WriteAlign(long long* edit, AlignmentSummary &al_sum, unsigned long Len)
{
  unsigned long i=0;
  unsigned long j=0;
  unsigned long k=0;
  unsigned long ac=0;
  unsigned long tmp;
  bool prev_del = false;
  mpq_class cost;
  while (i<Len)
  {
    if (edit[i]==1)
    {
      ASeqChar1[ac]=Seq1[j];
      ASeqChar2[ac]=Seq2[k];
      (*al_sum.matrix_summary)[make_pair(Seq1[j], Seq2[k])]+=1; 
      j++;
      k++;
      ac++;
      prev_del=false;
    }
    else if (edit[i] < -1)
    {
      tmp = k;
      al_sum.GapStart += 1;
      for (k; k < tmp-edit[i]-1; k++)
      {
        al_sum.Gap += 1;
        ASeqChar1[ac]='-';
        ASeqChar2[ac]=Seq2[k];      
        ac++;
      }
      al_sum.Gap -= 1;
      prev_del=false;
    }
    else if (edit[i] > 1)
    {
      tmp = j;
      if (!prev_del)
      {
        al_sum.GapStart += 1;
        al_sum.Gap -= 1;
      }
      for (j; j < tmp+edit[i]-1; j++)
      {
        al_sum.Gap += 1;
        ASeqChar1[ac]=Seq1[j];
        ASeqChar2[ac]='-';
        ac++;
      }
      prev_del=true;
    }
    i++;
  }
  ASeqChar1[ac]='\0';
  ASeqChar2[ac]='\0';
  printf("%s\n\n%s\n",ASeqChar1,ASeqChar2);
  cout << "Score: " <<  al_sum.Score() << endl;
  cout << al_sum.PrintSummary() << endl;

  return cost;
}


mpq_class WorkBench::SubMatrix(char i, char j)
{
  return subMatrix[param_idx[i]*param_idx.size()+param_idx[j]];
}

WorkBench::WorkBench(AjPSeq query, AjPSeq dbSeq)
{
  Sequence1 = query;
  Sequence2 = dbSeq;
  Seq1 = (char *) ajStrGetPtr(ajSeqGetSeqS(Sequence1));
  Seq2 = (char *) ajStrGetPtr(ajSeqGetSeqS(Sequence2));
  Seq1_len = ajSeqGetLen(Sequence1);
  Seq2_len = ajSeqGetLen(Sequence2);
  ASeqChar1 = new char[Seq1_len+Seq2_len];
  ASeqChar2 = new char[Seq1_len+Seq2_len];
  ASequence1 = ajSeqNewRangeC(ASeqChar1, ajSeqGetOffset(Sequence1), 
               ajSeqGetOffend(Sequence1), ajSeqIsReversed(Sequence1));
  ASequence2 = ajSeqNewRangeC(ASeqChar2, ajSeqGetOffset(Sequence2),
               ajSeqGetOffend(Sequence2), ajSeqIsReversed(Sequence2));

}

void WorkBench::PrintParams()
{
  cout << "Gap Start: " << gapStart << " Gap: " << gap << endl << endl;
  cout << "Substitution Matrix: " << endl << "  ";
  cout << endl;
  for(long i=0; i<param_idx.size(); i++)
  {
    for(long j=0; j<param_idx.size(); j++)
      cout << " " << subMatrix[i*param_idx.size()+j].get_str();
    cout << endl;
  }
}




AjPMatrixf WorkBench::FloatSubs()
{
  AjOMatrixf* matrix = 0;  
  //AJNEW0(matrix);
  matrix = new AjOMatrixf;
/*
AjPMatrixf    ajMatrixfNew (const AjPPStr codes, ajint n,
                            const AjPStr filename);
*/

//  char* acodes[] = {"A",   "T",   "G",   "C",   "S",   "W",   "R",   "Y",   "K",   "M",   "B",   "V",   "H",   "D",   "N",   "U"};
//  matrix = ajMatrixfNew((const AjPPStr) acodes, 16, (const AjPStr)"params.txt");
  matrix->Size=param_idx.size()*param_idx.size();
  matrix->SizeRow=param_idx.size();
//  #ifdef DEBUG3
//  cout << "matrix->SizeRow: " << matrix->SizeRow << endl; 
//  #endif
  char codes[matrix->SizeRow];
/*  matrix->Codes = new AjPStr[matrix->SizeRow];
  matrix->Codes[0]=(AjPStr) "ab";
  matrix->Codes[1]=(AjPStr) "ab";
  matrix->Codes[2]=(AjPStr) "ab";
  matrix->Codes[3]=(AjPStr) "ab";
  matrix->Codes[4]=(AjPStr) "ab";
  matrix->Codes[5]=(AjPStr) "ab";
  matrix->Codes[6]=(AjPStr) "ab";
  matrix->Codes[7]=(AjPStr) "ab";
  matrix->Codes[8]=(AjPStr) "ab";
  matrix->Codes[9]=(AjPStr) "ab";
  matrix->Codes[10]=(AjPStr) "ab";
  matrix->Codes[11]=(AjPStr) "ab";
  matrix->Codes[12]=(AjPStr) "ab";
  matrix->Codes[13]= (AjPStr) "ab";
  matrix->Codes[14]=(AjPStr) "ab";
  matrix->Codes[15]=(AjPStr) "ab";  */

  char c2[2];
  //float vals[matrix->SizeRow*matrix->SizeRow]; 
  AjFloatArray vals[matrix->SizeRow*matrix->SizeRow];
  map<char,int>::iterator mpr = param_idx.begin();
  map<char,int>::iterator mpc = param_idx.begin();


  //AJCNEW0(matrix->Codes, matrix->SizeRow);
  //int i;
  //for(i=0; i<matrix->SizeRow; i++)
  //    matrix->Codes[i] = ajStrNew();

  //AJCNEW0(matrix->Matrixf, matrix->SizeRow);
  //for(i=0; i<matrix->SizeRow; i++)
  //    AJCNEW0(matrix->Matrixf[i], matrix->SizeRow);
  //matrix->Cvt = ajSeqcvtNewStr(codes, matrix->SizeRow);
  map<char,int>::iterator test  = param_idx.end();
  --test;
  cout << "param start: " <<  (*mpr).second << endl;
  cout << "param end: " <<  (*test).second << endl;
  PrintParams(); 
  while (mpr != param_idx.end())
  {
    codes[(*mpr).second]=(*mpr).first;
    float val_row[matrix->SizeRow];
    while (mpc != param_idx.end())
    {
      cout << (*mpr).second << " " << matrix->SizeRow << " " << (*mpc).second << endl;
      //vals[(*mpr).second * matrix->SizeRow + (*mpc).second]=SubMatrix((*mpr).first, (*mpc).first).get_d();
      
      val_row[(*mpc).second]=SubMatrix((*mpr).first, (*mpc).first).get_d();
      ++mpc;
    }
    vals[(*mpr).second]=val_row;
    ++mpr;
  }
  matrix->Name = ajStrNewC("_____");
//  matrix->Matrixf=vals;
  for(int i=0; i<matrix->SizeRow; i++)
  {
    c2[0]=codes[i];
    c2[1]='\0';
    //matrix->Codes[i] = new AjOStr;
    ajStrAssignC(&matrix->Codes[i], c2);

  }
  return matrix;

}

void WorkBench::ReadParams(string file_name)
{
  ifstream in;
  in.open(file_name.c_str());
  char line[2000];
  string word; 
  while(!in.eof())
  {
    in >> word;
    if (word == "GAP")
    { 
      word = ""; in >> word;
      gapStart = word;
      word = ""; in >> word;
      gap = word;       
    }
    else if (word=="SUB")
    {
      in.getline(line, 2000);
      in.getline(line, 2000);
      int i=0;
      char *tok;
      tok = strtok(line," \t");
      while (tok != NULL)
      {
        param_idx.insert(make_pair(toupper(tok[0]),i)); 
        tok = strtok(NULL, " \t");
        i++;
      }
      subMatrix = new mpq_class[param_idx.size()*param_idx.size()];
      i=0;
      while(!in.eof())
      {
        in.getline(line, 2000);
        tok = strtok(line," \t");
        while(tok != NULL)
        {
          subMatrix[i]=tok; 
          tok = strtok(NULL, " \t");
          i++;
        }
      }
    }  
  }
  return;
}
