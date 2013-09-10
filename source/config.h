
//#define DEBUG_GLOBAL

/* DEBUG1 is used for:
   1)Gotoh phase debugging

*/
//#define DEBUG1 //Not thread safe

/*
  DEBUG2 is used for:
  1) Debugging the M&M algorithm,
  specifically MainLoop and Diff

  //Note: Requires DEBUG1 at this time.
*/

//#define DEBUG2

/*
DEBUG3 used for debugging ListWalk and
alignment string creation.
*/

//#define DEBUG_GLOBAL

#define OMP_ACTIVE

#ifdef DEBUG_GLOBAL
  #define DEBUG1
  #define DEBUG2
  #define DEBUG3
#endif
