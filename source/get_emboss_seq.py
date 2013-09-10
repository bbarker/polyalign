import re
import sys


def seqname(sn):
  """
  returns truncated sequence name.
  """
  return sn[0:6]

def getseq(tsn):
  seq = ""
  alignment = file(sys.argv[2],"r")
  for line in alignment.readlines():
    if line[0:6]==tsn:
      seq += line[7:57]
  alignment.close()
  return seq

tsn = seqname(sys.argv[1])
seq = getseq(tsn)
print seq
