polyalign
=========

Parallel and Parametric Alignment tools are primarily designed for
exploring questions about alignment methods, though they may also be
useful to biologists. https://code.google.com/p/polyalign/ is the
original svn repository (for historical purposes).

Emboss problems mail emboss-bug@embnet.org

1) Extract emboss archive, Go to the emboss directory 

cd emboss

(build emboss)

./autogen.sh

./configure

make

Note that even if the build fails prematurely, it may be
possible enough completed to proceed with subsequent steps
below.

2) make new directory embassy if it does not exist already.

mkdir embassy

3) Go into that directory

cd embassy

4) If you downloaded a zip file:

unzip Download_Location/polyalign-master.zip

mv polyalign-master polyalign

Otherwise, just clone the repository:

git clone https://github.com/bbarker/polyalign.git

5) go into the polyalign directory:

cd polyalign

6) configure and compile

./autogen.sh (this may be all that is needed)

./configure (use same options as you used to compile emboss)

make

7) To test:

cd source

para2 -asequence seqa.fasta -bsequence seqb.fasta

(press enter)


======================================================================
Make sure your compiler supports OpenMP.  GCC 4.2.0+ is an OpenMP
compiler that is freely available.

Set the following environmental variables prior to running the
application:

OMP_NUM_THREADS=Number of threads to use.

Typically this is the number of processing cores on your system.  Note
that this only works on shared memory systems and cannot benefit from
clusters and the like.

OMP_NESTED=TRUE

======================================================================

