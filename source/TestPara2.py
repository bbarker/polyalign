#!/usr/bin/python

import sys
import random
import math
import os
import time
import subprocess
import re

import pprint
pp = pprint.PrettyPrinter(indent=4)

def print_help():
  print "syntax is: randseq2.py {protX,prot,dna,dnan} length"

dna = dict({0:'A',1:'G',2:'C',3:'T'});
dnan = dict({0:'A',1:'G',2:'C',3:'T',4:'N'})
prot = dict({0:'G',1:'A',2:'L',3:'M',4:'F',5:'W',6:'K',7:'Q',8:'E',9:'S',10:'P',11:'V',12:'I',13:'C',14:'Y',15:'H',16:'R',17:'N',18:'D',19:'T'})
protX = dict({0:'G',1:'A',2:'L',3:'M',4:'F',5:'W',6:'K',7:'Q',8:'E',9:'S',10:'P',11:'V',12:'I',13:'C',14:'Y',15:'H',16:'R',17:'N',18:'D',19:'T',20:'X'})

def random_seq(seqtype, length):
  seq=">seq\n"
  if seqtype=="dna":
    ldict=dna;
  elif seqtype=="dnan":
    ldict=dnan;
  elif seqtype=="prot":
    ldict=prot;
  elif seqtype=="protX":
    ldict=protX;
  else:
    print_help()
    raise SystemExit
  for i in range (0,length):
    seq += ldict[random.randint(0,len(ldict)-1)]
  seq+="\n"
  return seq

def random_jc_params():
  match=random.randint(5,10)
  mismatch=match-random.randint(0,5)
  gapstart=match-random.randint(5,15)
  gapextend=match-random.randint(5,15)
  return (match,mismatch,gapstart,gapextend)

def write_jc_matrices():
  #first for Para2
  jc_params = random_jc_params();
  para2_text = "GAP\t"+str(jc_params[2])+"\t"+str(jc_params[3])+"\n\nSUB\n"
  para2_text += "A\tC\tT\tG\n"
  for i in range (0,4):
    for j in range (0,4):
      if i==j:  para2_text += (str(jc_params[0])+"\t")
      else:     para2_text += (str(jc_params[1])+"\t")
    para2_text += "\n"
  para2_mat_file = open("params.txt","w")
  para2_mat_file.write(para2_text)
  para2_mat_file.close()
  #now for EMBOSS
  emboss_text="    A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U\n"
  emboss_chars=('A','T','G','C','S','W','R','Y','K','M','B','V','H','D','N','U')
  for i in range (0,16):
    emboss_text+=(emboss_chars[i]+"   ")
    for j in range (0,16):
      if i==j:  emboss_text += (str(jc_params[0])+"   ")
      else:     emboss_text += (str(jc_params[1])+"   ")
    emboss_text+="\n"
  #print emboss_text
  emboss_mat_file = open("DNAMM","w")
  emboss_mat_file.write(emboss_text)
  emboss_mat_file.close()
  #print para2_text
  #print emboss_text
 
  return (jc_params[2],jc_params[3])

log_mismatch = open("log.txt","w")
log_mismatch.write("Max Seq Length: " + sys.argv[1] + "\n")
while(1):
  seqA = random_seq("dna",int(random.randint(1,int(sys.argv[1]))))
  seqf = open("seqa.fasta","w")
  seqf.write(seqA)
  seqf.close()
  seqB = random_seq("dna",int(sys.argv[1]))
  seqf = open("seqb.fasta","w")
  seqf.write(seqB)
  seqf.close()
#turn off random matrices now:
  for i in range (1,500):
    log_mismatch.flush()
    gaps = write_jc_matrices()
    stretcher_cmd="../../../emboss/stretcher -datafile DNAMM -gapopen "+str(-gaps[0])+" -gapextend "+str(-gaps[1])+" -asequence seqa.fasta -bsequence seqb.fasta"
    #stretcher_cmd="/home/brandon/emboss/emboss/emboss/stretcher -datafile DNAMM -gapopen 4 -gapextend 1 -asequence seqa.fasta -bsequence seqb.fasta"
    #ZZ(stretcher.stdin,stretcher.stdout) = os.popen4(stretch_cmd,"w")
    stretcher = subprocess.Popen(stretcher_cmd, shell=True, bufsize=1000000, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    stretcher.stdin.write("\n\n")
    stretcher.stdin.flush()
    stretcher.wait()
    para2_cmd="../source/para2.go -asequence seqa.fasta -bsequence seqb.fasta"
    #(para2_in,para2_out) = os.popen4(para2_cmd,"w")
    para2 = subprocess.Popen(para2_cmd, shell=True, env=para2_env, bufsize=1000000, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    para2.stdin.write("\n")
    para2.stdin.flush()
    para2.wait()
    #print para2_cmd
    #print stretcher_cmd
    #read in output
    #stretcher first:
    stretcher_file_out = open("seq.stretcher","r")
    stretcher_output = stretcher_file_out.read()
    #para2:
    para2_output=para2.stdout.read(-1)
    #Rvar_rex = re.compile('([a-zA-Z_]+\w*)($|\+|\-|\s|\*|\/|\^|(|))')
    #p_parse = Rvar_rex.findall(repr(p))
    #Rvars = dict([(var,'0') for (var,d1,d2) in p_parse])
    #get scores
    score_rex = re.compile('(Score:\s+)(\-?\w+)')
    #stretcher:
    stretcher_score = score_rex.findall(stretcher_output)
    #para2:
    para2_score = score_rex.findall(para2_output)
    #pp.pprint(stretcher_score)
    #print "Score: " + stretcher_score[0][1] + " : " + para2_score[0][1]
    if int(stretcher_score[0][1]) != int(para2_score[0][1]):
      print "Score: " + stretcher_score[0][1] + " : " + para2_score[0][1]
      log_mismatch.write(seqA)
      log_mismatch.write(seqB)
      log_mismatch.write(stretcher_cmd)
      log_mismatch.write(para2_cmd)
      log_mismatch.write(stretcher_output)
      log_mismatch.write(para2_output)
      dnamm_file = open("DNAMM","r")
      dnamm = dnamm_file.read()
      log_mismatch.write(dnamm)
      log_mismatch.write("######################################################################\n\n\n\n\n")
      log_mismatch.flush()

    stretcher.stdin.close()
    stretcher.stdout.close()
    stretcher_file_out.close()
    para2.stdin.close()
    para2.stdout.close()
  
  #cleanup






