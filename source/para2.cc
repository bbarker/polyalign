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
#include <iostream>
#include "alignment.h"

using namespace Para2;
int main(char argc, char **argv)
{
  AjPSeqall seqall;
  AjPSeq SeqA;
  AjPSeq SeqB;
  AjPAlign align;
  ajint nd;

  embInit ("para2", argc, argv);

  SeqA = ajAcdGetSeq("asequence"); 
  ajSeqTrim(SeqA);
  ajSeqFmtUpper(SeqA);
  seqall = ajAcdGetSeqall("bsequence");
  ajSeqallNext(seqall,&SeqB);
  ajSeqTrim(SeqB);
  ajSeqFmtUpper(SeqB);
  align = ajAcdGetAlign("outfile");
  GlobalAlign test(SeqA, SeqB);
  test.MainLoop();
  //test.PrintParams();
  //test.score = test.ListWalk(test.optimal_vertices, *test.alignment_summary);
  test.alignment_summary->PrintSummary();
  test.score = test.WriteAlign(test.Edits, *test.alignment_summary, test.Seq1_len+test.Seq2_len);
  test.alignment_summary->PrintSummary();
  //cout << test.score << endl;

//Need to get EMBOSS integration fully working.

/*  ajAlignDefineSS(align, test.ASequence1, test.ASequence2);
  ajAlignSetGapR(align, test.gap.get_d(), test.gapStart.get_d());
  //This of course may not display correctly.
  ajAlignSetScoreR(align, test.score.get_d());
  ajAlignSetMatrixFloat(align, test.FloatSubs());
  ajAlignSetStats(align, -1, strlen(test.ASeqChar1), nd, -1, -1, NULL);

  ajAlignWrite(align);
  ajAlignClose(align);
*/
  ajAlignDel(&align);

  ajSeqDel(&SeqA);
  ajSeqDel(&SeqB);
  ajSeqallDel(&seqall);

  embExit();
  return 0;
}

