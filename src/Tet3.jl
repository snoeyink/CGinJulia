""" Tet3
Tetrahedra are stored as lex min order that have positive orientation.
"""
module Tet3
jd,Coord} # w,x,y,z,1
struct Sphere3
	w::Coord4
	x::Coord4
	y::Coord4
	z::Coord4
end
function Sphere3(p::Point3L,q::Point3L,r::Point3L,s::Point3L)
	qp = q.-p
	rp = r.-p
	sp = s.-p
	Sphere3()
end

dot(p::Point3L, s::Sphere3) = s.w+s.x*p.x+s.y*p.y+s.z*p.z+p.q



## incremental delaunay statuse structure
const TETS_PER_VERTEX = 12
const CORNERS_PER_VERTEX = 4TETS_PER_VERTEX
struct Del3d
	vert::Vector{Point3L} # list of input points; first must be at infinity and first 5 are affine independent
	sph::Vector{Sphere3} # sphere for each tetrahedron
	s::Vector{Corner} # corner table: 4 per tetrahedron with vertices ordered lex min positive
	cornerlists::MVector{Corner[0,-1]} # head of free list free and live corner indices
	function Del3d(vert)
		vert[1] != Inf3L && pushfirst!(vert, Inf3L)
		N = length(vert)
		sph = Vector{}
		new(vert, sph, [Corner(i+3) for i = 1:CORNERS_PER_VERTEX*N], @MVector Corner[0,-1] )
end

TETRA(corner) = (corner >> 2)
CORNER(tet, index) = (corner<<2 + index)
BASECORNER(corner) = (corner & 0xFFFFFFFC) #zero last two bits
LASTCORNER(corner) = (corner | 3)

## making vertices

rand(0:0x3fff, 3)

## stacks
# dfs: all dead
dfs = Corner[] # stack for depth first earch of neighboring tetrahedra
idfs = Corner[] # stack for tetras adjacent to infinite vertex
nhbr = Corner[] # stack for dead corners with live neighbors
kill = Corner[] # stack for base corners of tetrahedra to kill

#= stack operations
push!(dfs,2,3)
length(dfs)
pop!(dfs)
resize!(dfs, 0)
length(dfs)
=#



## Tetrahedron manipulation tables
# c+offset[i,1+(c mod 4)] advances c to (c+i)mod5
#INDEX(c::Int) = mod(c-1,4)+1
const mask3 = Int32(3)
@inline INDEX(c::Int) = ((c-1)&mask3)+1
const offset = @SMatrix  Int32[0 0 0 0; 1 1 1 -3; 2 2 -2 -2; 3 -1 -1 -1]
@inline INCREMENT(c::Int) = @inbounds (c+offset[1,1 + (c & mask3)])
@inline INC(c::Int) = @inbounds (c+offset[1,mod4a(c)])
##
#drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side.
const drop =  @SMatrix Int32[ 3  2  4;  1  3  4;  2  1  4;  1  2  3]
# offdr[i] contains drop(i)-index(i)
const offdr = @SMatrix Int32[  2   1   3;  -1   1   2;  -1  -2   1;  -3  -2  -1]
# invdrop[i][k] = j whenever drop[i][j] = k.  4s signal i=k; bad because i is dropped.
const invdrop = @SMatrix Int32[ 5  2  1  3;  1  5  2  3;  2  1  5  3;  1  2  3  5]

function checkInc(N::Int)
  y = zero(Int32)
  for i in 1:N
       y += INC(i)
  end
end
# @btime Tet3.checkInc(100_000) with these different mods in INC
# I'm surprised that a conditional is better.
mod4a(c::Int) = mod1(c,Int32(4)) #  125.299 μs (0 allocations: 0 bytes)
mod4b(c::Int)::Int32 = ((c-one(typeof(c)))&mask3)+one(Int32) # slower:   174.299 μs (0 allocations: 0 bytes)
mod4c(c::Int)::Int32 = mod(c-one(typeof(c)),Int32(4))+one(Int32) # slower:  174.299 μs (0 allocations: 0 bytes)
#@inline INC(c::Int) = @inbounds (c+offset[one(Int32),one(Int32)+(c & mask3)]) #   125.299 μs (0 allocations: 0 bytes)


## d3permute.h
const MASKTAIL = 0x3F # mask for bits in tail (for levels)
const NUMTAIL  = 0x40 # how many bit patterns in tail
const NUMTAIL2 = 0x80 # 2*NUMTAIL

# coordinate info to record while reading
bboxType = Int64[6];

# Bounding box & coord histogram functions
# WARNING: I'm being lazy here, and using globals
extern bboxType bb;  #- global: assume that we open one file at a time. -#
const BBMAXCOORD 0xFFFffff # 28 bits should be enough
const bbInit() { bb[0] = bb[1] = bb[2] =  BBMAXCOORD; #-min x,y,z-#\
                   bb[3] = bb[4] = bb[5] = -BBMAXCOORD; #-max x,y,z-#}

# -DNOCOORDHIST will turn of the use of coordinate tail histograms for choosing levels.
#ifdef NOCOORDHIST
const cHistInit()
const LEVELNO(pv) levelLUT[((((int)(pv->x)) | ((int)(pv->y)) \
 	            | ((int)(pv->z)) ) & MASKTAIL)]
#else
typedef short coordHistType[3][NUMTAIL2]; # histograms for last 7 coord bits
extern coordHistType cHist; #- global: assume that we open one file at a time. -#
const cHistInit() { bzero((void *)cHist, sizeof(cHist)); }
const LEVELNO(pv) levelLUT[(((xorBits[XX]^((int)(pv->x))) | (xorBits[YY]^((int)(pv->y))) \
 	            | (xorBits[ZZ]^((int)(pv->z))) ) & MASKTAIL)]
#endif

# Called by readers to set coordinates. WARNING: uses global variables bb and cHist
const TWO26 0x2000000 # 2^26 added to coords to make them positive
inline void setVert(ppointType v, int indx, double xx, double yy, double zz,
	     double rad, double mult) {  # this part of the code may be I/O bound anyway,
  double mr = mult*rad;                  # so we count most common coordinate bit tails
  int ix, iy, iz, i;
  (v)->index = indx;
  (v)->x = ix = (int)(mult*(xx)+0.5);
  (v)->y = iy = (int)(mult*(yy)+0.5);
  (v)->z = iz = (int)(mult*(zz)+0.5);
  (v)->sq = (int)(-mr*mr);

  if (bb[0] > ix) bb[0] = ix; # update bounding box
  if (bb[1] > iy) bb[1] = iy;
  if (bb[2] > iz) bb[2] = iz;
  if (bb[3] < ix) bb[3] = ix;
  if (bb[4] < iy) bb[4] = iy;
  if (bb[5] < iz) bb[5] = iz;

#ifndef NOCOORDHIST
  i = (ix+TWO26)&MASKTAIL; cHist[XX][i]++; # add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; # 5 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; # 4 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; # 3 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; # 2 bits
  i = (i>>1)+NUMTAIL; cHist[XX][i]++; # 1 bit

  i = (iy+TWO26)&MASKTAIL; cHist[YY][i]++; # add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; # 5 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; # 4 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; # 3 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; # 2 bits
  i = (i>>1)+NUMTAIL; cHist[YY][i]++; # 1 bit

  i = (iz+TWO26)&MASKTAIL; cHist[ZZ][i]++; # add to coordinate histograms: 6 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; # 5 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; # 4 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; # 3 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; # 2 bits
  i = (i>>1)+NUMTAIL; cHist[ZZ][i]++; # 1 bit
#endif
}


# permutation by boxOrder
# box indices
const BX(x,mask) (((x)>>bitshift)&mask)
typedef int boxMatrixType[8][8][8]; # a boxMatrix has an int for each 3 bits of x,y,z
typedef int *pboxMatrixEntryType; # pointer to a box matrix entry

#increment boxMatrix entry (based on bitshift & mask)
const boxv(bm, v, mask) (bm[BX((int)((v)->x),mask)][BX((int)((v)->y),mask)][BX((int)((v)->z),mask)])#--#
const boxIncr(bm, x,y,z, mask) (bm[BX(x, mask)][BX(y, mask)][BX(z, mask)]++)
const boxIncrv(bm, v, mask) boxIncr(bm, (int)((v)->x), (int)((v)->y), (int)((v)->z), mask)

## d3permute.c

rr#- d3permute.c                Jack Snoeyink July 2003
Using the various readers (pdbreader, plyreader, etc),
we read in coordinates and compute bounding boxes
and levels.  -#

#include <time.h>
#include "d3permute.h"

# boxOrder gives the order of boxMatrix elements for the Hilbert curve generators associated
# with the 12 different edges.  hilbertCase gives the generator number for boxOrder[0].
#include "d3boxOrder.c"

#I'm being lazy here and using some globals for coordinate information
bboxType bb;  #- global: assume that we open one file at a time. -#
#ifndef NOCOORDHIST
coordHistType cHist; #- global: -#
#endif

boxMatrixType bm, bm2;  #- global: two box levels max -#


const NLEVELS 7
# NLEVELS-1 - # of trailing zeros is level number
const levelLUT = @SVector Int8[0 6 5 6 4 6 5 6 3 6 5 6 4 6 5 6 2 6 5 6 4 6 5 6 3 6 5 6 4 6 5 6 1 6 5 6 4 6 5 6 3 6 5 6 4 6 5 6 2 6 5 6 4 6 5 6 3 6 5 6 4 6 5 6]
# for BBbits
const high = @SVector Int8[0  1  2 2  3 3 3 3  4 4 4 4 4 4 4 4  5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5]

inline int BBbits(bboxType bb) { # returns the number of bits in bb range
  int range;
  int mask = 0xffffFFE0; # mask last 5 bits
  int shift = 0;
  range = bb[1]-bb[0]; # find max range in all three coords
  if (range < bb[3]-bb[2]) range = bb[3]-bb[2];
  if (range < bb[5]-bb[4]) range = bb[5]-bb[4];

  while (range & mask) { # some ones bits left above last 5.
    mask <<= 5;
    shift += 5;
  }
  return shift+high[range>>shift]; # shift needed to make all bits zero
}

# This function permutes points from v to newv by boxOrder
# From a coordinate x, it uses bits (x>>bitshift)&mask.  nboxes = (mask+1)^3;
# hcase is so we can handle hilbert cases for first level.
# This function is here for documentation.  We don't actually call it,
#  but weave this code into the reading, tuning it for efficiency.
#
#-
inline void permute(boxMatrixType bm, const ppointType v, const int nvert, ppointType newv,
		    const short bitshift, const short mask, const short nboxes, const short hcase) {
  pboxMatrixEntryType pbm;
  ppointType pv;
  int i, n, tmp;

  bzero((void *)(bm), sizeof(bm)); # init box counts
  for (pv = v; pv < v+nvert; pv++)
    boxIncrv(bm, pv, mask); # count for each box

  n = 0; # prefix sum the boxes in boxOrder
  for (i = 0; i < nboxes; i++) {
    pbm = ((pboxMatrixEntryType) bm) + boxOrder[hcase][i];
    tmp = *pbm;
    *pbm = n;
    n += tmp;
  }

  ASSERT(n == nvert, "total wrong after bm assignment");
  for (pv = v; pv < v+nvert; pv++)
    newv[boxIncrv(bm, pv, mask)] = *pv; # move to place
}
#--#
void bbPrint() {
  int i;
  printf("BB %d (%f %f %f; %f %f %f)\n", BBbits(bb),
	 bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
#ifndef NOCOORDHIST
  for (i = 0; i<NUMTAIL2; i++)
    printf("%5d %5o  %4d %4d %4d\n", i, i, cHist[XX][i], cHist[YY][i], cHist[ZZ][i]);
#endif
}

const NBOXES 64
# read vertices & return coordinfo
void d3permute(FILE *fid, int (*reader)(FILE *, ppointType),
	    ppointType *vout, int *nvert) {
  int i, j, k, tmp;
  int bitshift, hcase;
  int n, nv, nvi;
  int nblock, sblock, eblock;
  int origx, origy, origz;
#ifndef NOCOORDHIST
  int xorBits[3];
#endif
  clock_t tic, toc; #--#
  ppointType vin, v;
  ppointType pv, pvi;
  pboxMatrixEntryType pbm;
  int nlev[NLEVELS]; # Number/offset per level

  bbInit(); # initialize bounding box
  cHistInit(); # initialize coordinate histogram
  vin = (ppointType) calloc(MAXVERT, sizeof(pointType)); # pdb has a five digit atom# field
  if (vin == NULL) { printf("ERROR: could not allocate memory for MAXVERT points\n"); exit(EXIT_FAILURE); }

  nvi = reader(fid, vin); # read in vertices & compute bbox and cHist

  tic = clock();#--#
  if (nvi > MAXVERT) { printf("FATAL ERROR: Too many points %d; increase MAXVERT\n", nvi); exit(EXIT_FAILURE); }
  vin = (ppointType) realloc(vin, (nvi+1)*sizeof(pointType));
  v = (ppointType) calloc(nvi+1, sizeof(pointType)); # one more for infinite point
  if (v == NULL) { printf("ERROR: could not allocate memory for %d points\n", nvi+1); exit(EXIT_FAILURE); }
  v[0].index = -1; v[0].x = v[0].y = v[0].z = 0; v[0].sq = 1;
  *nvert = nvi+1; # output results
  *vout = v;

#ifndef NOCOORDHIST
  for (k = XX; k<= ZZ; k++) { # i,i+1 are bit patterns we decide between
    i = (MASKTAIL-1)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 1 bit
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 2 bit
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 3 bit
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 4 bit
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 5 bit
    i = (i&MASKTAIL)<<1; if (cHist[k][i] < cHist[k][i+1]) i++; # 6 bit
    xorBits[k] = i; # frequent bit pattern for this coordinate
  }
#endif

  bzero((void *)(nlev), sizeof(nlev)); # zero, then accumulate level counts for vin
  for (pvi = vin; pvi < vin+nvi; pvi++)
    nlev[LEVELNO(pvi)]++;

  n = 1;
  for (i = 0; i<NLEVELS; i++) { # prefix sum level counts for level offsets
    tmp = nlev[i];
    nlev[i] = n;
    n += tmp;
  }
  ASSERT(n==1+nvi, "bad count after level number prefix sum");

  for (pvi = vin; pvi < vin+nvi; pvi++) { # copy levels to v; adjust coords to 0--range
    pv = v + (nlev[LEVELNO(pvi)]++);
    VSUBASSN(pvi, bb[0], bb[0], bb[0], pv);
  }

  bitshift = BBbits(bb)-3; # find how many bits we need for entire range

  nv = 1; # index of next point to consider in v
   printf("Levels ");
  for (j = 0; j < NLEVELS; j++) { # order each level
    #-    printf("Level %d has %d [%d,%d): ", j, nlev[j]-nv, nv,nlev[j]);  #--#
    printf("%d:%d ", j, nlev[j]-nv);  #--#
    if (nlev[j]-nv > 64) { # 512 boxes on this level:  v[nv..nlev[j]) permuted using vin[nv..nlev[j])

      bzero((void *)(bm), sizeof(bm)); #- init box counts -#
      for (pv = v+nv; pv < v+nlev[j]; pv++)
	boxIncrv(bm, pv, 7); # count vertices of v[nv..nlev[j]) in each box
      n = nv;

      for (i = 0; i < 512; i++) { #- prefix sum the boxes in boxOrder -#
	#-printf("<%d:%o ", n, i);#--#
	pbm = ((pboxMatrixEntryType) bm) + boxOrder[j&1][i]; # alternating directions on levels
	tmp = *pbm; *pbm = n; n += tmp;
      }
      #-printf("<%d]\n", n);#--#

      ASSERT(n == nlev[j], "total wrong after bm assignment");
      #-      printf("\n L%d:"); #--#
      for (pv = v+nv; pv < v+nlev[j]; pv++) {
	#-	printf(" %d>%d,", pv->index, boxv(bm, pv, 7)); #--#
	vin[boxIncrv(bm, pv, 7)] = *pv; #  move to place in vin
      }
      # we now have this level in vin[nv..nlev[j]), and need it back in v[nv..nlev[j])
      # groups with few points can just be copied; those with many get radix sorted.
      # bm[] contains the prefix sum of counts INCLUDING current. I.e. old bm[1] is now bm[0]

      # We now consider radix sorting blocks in this level.
      # invariants: sblock = start of curr block; nvi = start of block that needs to be copied to v
      # i is block number, eblock is end of current block, nblock is # in block
      eblock = nvi = nv;
      for (i = 0; i < 512; i++) { # undo prefix sum (shifted) to get number in block
	sblock = eblock;
	eblock = ((pboxMatrixEntryType) bm)[boxOrder[j&1][i]];
	nblock = eblock - sblock;
	if (nblock > 32){ # shuffle this block
	  if (nvi < sblock) # There are old blocks to copy first
	    memcpy((void *)(v+nvi), (const void *)(vin+nvi), (sblock-nvi)*sizeof(pointType));

	  #-      printf("\n   c(%d,%d)",nvi, sblock); #--#
	  bitshift -= 3;
	  hcase = hilbertCase[j&1][i];
	  # Now we shuffle and copy vin[sblock..eblock) to v[sblock..eblock), with permutation

	  bzero((void *)(bm2), sizeof(bm2)); #- init box counts -#
	  for (pvi = vin+sblock; pvi < vin+eblock; pvi++)
	    boxIncrv(bm2, pvi, 7); # count vertices in each box
	  n = sblock;
	  for (k = 0; k < 512; k++) { #- prefix sum the boxes in boxOrder -#
	    #-printf("{%d:%o|%o ", n, i, k);#--#
	    pbm = ((pboxMatrixEntryType) bm2) + boxOrder[hcase][k];
	    tmp = *pbm;
	    *pbm = n;
	    n += tmp;
	  }
	  #-printf("{%d}\n", n);#--#
	  ASSERT(n == eblock, "total wrong after bm assignment");
	  for (pvi = vin+sblock; pvi < vin+eblock; pvi++) {
	    #-	    printf(" %d>%d,", pvi->index, boxv(bm2, pvi, 7)); #--#
	    k = boxIncrv(bm2, pvi, 7);
	    v[k] = *pvi;
	  }

	  bitshift += 3; # undo bitshift change above
	  nvi = eblock; # Here's where we be after shuffle & copy
	}
      }
      if (nvi < nlev[j]) # There are still old blocks to copy to finish this level
	memcpy((void *)(v+nvi), (const void *)(vin+nvi), (nlev[j]-nvi)*sizeof(pointType));
    }
    nv = nlev[j];
  }
  fflush(stdout);
  free(vin);

  origx = (bb[3]+bb[0])/2; # center the bbox
  origy = (bb[4]+bb[1])/2;
  origz = (bb[5]+bb[2])/2;
  # get rid of some duplicates as we lift
  pvi = v+1; # pointer to previous
  for (pv = pvi; pv < v+(*nvert); pv++) {
# move origin and lift (assumes point has been lifted by -rad^2)
     pv->x -= origx; pv->y -= origy; pv->z -= origz; pv->sq += DOT(pv,pv);
     if (!EQUALPV(pv,pvi)) # new point
       pvi++;
     else if (pv->sq >= pvi->sq) # keep old point
       continue;
     *pvi = *pv;
  }

  toc = clock();
  printf("Hilbert (%d -> %d pts) time(secs) %f\n", *nvert, pvi-v,
	 ((double) (toc - tic)) / CLOCKS_PER_SEC);  (void)fflush(stdout); #--#

#ifndef NOASSERT
  {
    short *index = (short *) calloc(*nvert, sizeof(short));
    bzero(index, (*nvert)*sizeof(short));
    for (i = 1; i < (pvi-v); i++)
      if (index[v[i].index] != 0)
	printf("POST: INDEX %d appears at %d and %d\n", v[i].index, i, index[v[i].index]);
      else
	index[v[i].index] = i;
    free(index);
  }
#endif

  *nvert = pvi-v;


}

## d3audit.c

# AUDIT routines for d3.c

#- Tetrahedra are groups of four corners in order of increasing
vertex index, except that the first two may be swapped to ensure
that the orientation determinant is positive.
I.e., take the lex smallest alternating permutation with positive sign.
There must always be an odd number of swaps between two permutations;
we swap the first two if necessary to achieve this.

The table indoff has the index offset for where
the vertex of index i will be found in the tetrahedron opposite c.
Vertex s[BASECORNER(c)+i].v will be at s[s[c].opp + indoff[INDEX(c)][INDEX(s[c].opp)][i]].v
(except that s[c].v and s[s[c].opp]].v are different, and opposite sides of common pl).
Note that indoff[*][j] uses each offset -j:4-j exactly once,
that indoff[i][j][i] = 0 for all possible i and j (so s[c].v is opposite of s[s[c].opp].v),
and that the tetra will never change, TETRA(s[c].opp) == TETRA(CORNERINOPP(i,c)).
-#

# This is the table of where the indices go; used in ASSERTS only.
# Replacing  \ With V falling at position c.opp:
#   corner c  \0ABCD 1ABCD 2ABCD 3ABCD
#       A    0: xxxx  BVCD  CBVD  BCDV
#       B    1: VACD  xxxx  ACVD  CADV
#       C    2: vBAD  AVBD  BAVD  ABDV
#       D    3: Vabc  bvac  ABVC  BACV
#
# Offsets from c.opp  I=-1, Z=-2, B = -3, H = -4
# Replacing \0ABCD 1ABCD 2ACBD 3ABCD
#     A    0: xxxx  0I12  0IZ1  0BZI
#     B    1: 1023  xxxx  Z0I1  Z0BI
#     C    2: 1203  I102  IZ01  BZ0I
#     D    3: 1230  1I20  ZI10  ZBI0
# Get these by subtracting 01234 from columns in order
const CORNERINOPP(i,c) (this->s[c].opp + indoff[INDEX(c)][INDEX(this->s[c].opp)][i])
const short indoff[4][4][4]
= {#-        0ABCD       1ABCD         2ACBD          3ABCD     -#
  #- 0: -# {{5,5,5,5},  { 0,-1, 1, 2},  { 0,-1,-2, 1},  { 0,-3,-2,-1}},
  #- 1: -# {{1,0,2,3},  { 5, 5, 5, 5},  {-2, 0,-1, 1},  {-2, 0,-3,-1}},
  #- 2: -# {{2,1,0,3},  {-1, 1, 0, 2},  {-1,-2, 0, 1},  {-3,-2, 0,-1}},
  #- 3: -# {{1,2,3,0},  { 1,-1, 2, 0},  {-2,-1, 1, 0},  {-2,-3,-1, 0}}};


void cornerPrint(const pd3stateType this, int c) {
  printf("%%%3d(%2d,%1d)%2d=", c, TETRA(c), INDEX(c), this->s[c].v-this->vert);
  if ((this->s[c].v >= this->vert) && (this->s[c].v < this->vert+MAXVERT))
    printf("(%d %5.0f %5.0f %5.0f) opp:%4d(%3d,%2d) \n", (infiniteV(this->s[c].v,this->vert)?0:1),
	   this->s[c].v->x, this->s[c].v->y, this->s[c].v->z, this->s[c].opp,
	   TETRA(this->s[c].opp),INDEX(this->s[c].opp));
  (void)fflush(stdout);
}

void cornerPrint4(const pd3stateType this, int c) {
  int b,j,k;
  psphereType sp;
  b = BASECORNER(c);
  if (this->sph != NULL) {
    sp = this->sph+TETRA(b);
    printf("disp('Sphere(%d) = <%5.0f %5.0f %5.0f %5.0f>')\n", TETRA(b),
	   sp->x,sp->y,sp->z,sp->sq);
  }
  printf("DetCheckH([");
  for (k = 0; k < 4; k++)
    printf(" %d %5.0f %5.0f %5.0f %5.0f; %% %d\n",
	   (infiniteV(this->s[b+k].v,this->vert)?0:1), this->s[b+k].v->x, this->s[b+k].v->y, this->s[b+k].v->z, this->s[b+k].v->sq, this->s[b+k].v->index);
  printf("]);\n");
  for (j=0;j<4;j++) cornerPrint(this, b+j);
}

# when keeping statistics...
int dropCnt =0;
int maxLocate =0, locateSideCnt = 0, sphereCnt = 0, startTetraCnt=0, freeTetraCnt=0, inSphereCnt=0, randbitCount=0;

#ifndef NOASSERT
const auditCorners(this, c) auditCornersAux(this, c);
#else #- NOASSERT -#
const auditCorners(this, c)
#endif #- NOASSERT -#

void auditCornersAux(const pd3stateType this, int SphereCheck) {
  int p, b,c, i,j,k, guard;
  double d;
  psphereType sp;
  ppointType vv, vp;

  for (p = 0; p < this->maxTetra; p++) {
    if (DEAD(p)) continue; # don't audit tetras on free list
    b = CORNER(p,0);
    for (c=b; c < CORNER(p,4); c++) { # per corner checks
      i = this->s[c].opp; # check opposite
      if (this->s[i].opp != c) { printf("%%AUDIT: wrong opp.opp \n"); cornerPrint(this, c); cornerPrint(this, i); }
      if (this->s[c].v == this->s[i].v){
	printf("%%AUDIT: Same vertex  s[%d(%d,%d)].v %d == opp[%d(%d,%d)].v %d \n",
	       c,TETRA(c), INDEX(c), this->s[c].v-this->vert, i, TETRA(i), INDEX(i), this->s[i].v-this->vert);
	cornerPrint4(this, c); cornerPrint4(this, i);
      }

      for (j = 0; j < 4; j++)
	if ((j != INDEX(c)) && (CORNERINOPP(j,c)<3))
	  if (this->s[BASECORNER(c)+j].v != this->s[CORNERINOPP(j,c)].v) {
	    printf("%%AUDIT:Bad auditCornerInOpp(%d,%d) = %d  since vertex %d != %d\n",
		   j,c,CORNERINOPP(j,c), this->s[BASECORNER(c)+j].v-this->vert, this->s[CORNERINOPP(j,c)].v-this->vert);
	    cornerPrint4(this, BASECORNER(c)); cornerPrint4(this, BASECORNER(CORNERINOPP(j,c)));
	    break;
	  }

      for (j = 0; j<4; j++) {
	k = CORNERINOPP(j,c);
	if ((TETRA(k) != TETRA(this->s[c].opp)) || (this->s[b+j].v != this->s[k].v)) {
	  if (TETRA(k) != TETRA(this->s[c].opp)) {
	    printf("%%AUDIT: CORNERINOPP(%d,%d) ==>%d: Accessing [%d:%d][%d:%d][%d]\n",
		   j, c, k-this->s[c].opp, c, INDEX(c), this->s[c].opp, INDEX(this->s[c].opp), j);
	    k = this->s[c].opp;
	  }
	  else {
	    if (b+j == c)
	      if (k == this->s[c].opp) continue; # these vertices are supposed to differ; don't flag them
	      else printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) shouldn't happen\n",
			  j, c, b+j, TETRA(b), j, k, TETRA(k), INDEX(k));
	    else {
	      printf("%%AUDIT: CORNERINOPP(%d,%d) says %d(%d,%d) and %d(%d,%d) should agree\n",
		     j, c, b+j, TETRA(b), j, k, TETRA(k), INDEX(k));
	      cornerPrint4(this, c);
	      cornerPrint4(this, b+j);
	      cornerPrint4(this, k);
	      ASSERT(TETRA(this->s[c].opp) == TETRA(CORNERINOPP(i,c)), "CORNERINOP screws up tetras");

	    }
	  }
	  cornerPrint4(this, b); cornerPrint4(this, BASECORNER(k)); break;
	}
      }

      # check sphere opposite corner c
      if (SphereCheck>0) 	# check sphere opposite corner
	{
	  k = this->s[c].opp;
	  sp = this->sph+TETRA(k);
	  if (infiniteV(this->s[c].v,this->vert)) {
	    d = spdotInf(sp, this->s[c].v);
	  } else {
	    vp = this->s[CORNER(TETRA(k),3)].v; # subtract from this
	    d = spdot(sp, this->s[c].v, vp);
	  }
	  if (d < 0) {
	    printf("disp('AUDIT: corner %d v%d in sphere %d(%d) =%5.0f');\n",
		   c, this->s[c].v-this->vert, TETRA(k), k, d);
	    cornerPrint(this, c); # print vertex
	    cornerPrint4(this, k);
	  }
	}
    }

    if (SphereCheck>0) 	# check sphere sqs (orient dets)
      {
	sp = this->sph+p; # only spheres using pt at infty have sq==0; none have sq < 0.
	b = CORNER(p,0);
	if (sp->sq < 0 || (sp->sq ==0 && !infiniteV(this->s[b].v,this->vert) && !infiniteV(this->s[b+1].v,this->vert))) {
	  printf("disp('AUDIT: sq<=0 in tetra %d(%d) =%5.0f'); \n", CORNER(p,0), p, sp->sq);
	  cornerPrint4(this, b);
	} #--#
	#-	if (pv-this->vert > 860 || ((pv-this->vert) % 100 == 0))
		for (vv = this->vert; vv < pv; vv++) { # Check all vertices against all spheres
		if (vv==this->s[b+0].v) continue;
		if (vv==this->s[b+1].v) continue;
		if (vv==this->s[b+2].v) continue;
		if (vv==this->s[b+3].v) continue;
		if (vv==this->s[b+4].v) continue;
		d = spdot(sp, vv);#	  d = spdot(sp, this->s[c].v);
		if (d < 0) {
		printf("disp('AUDIT: vertex v%d in sphere %d(%d) =%5.0f');\n", vv-this->vert, p, b, d);
		printf("   (%5d; %5.0f %5.0f %5.0f; %10f)\n",
		vv->index, vv->x, vv->y, vv->z, vv->sq);
		cornerPrint4(this, b);
		}
		}#--#
      }
  }
}

## bits.c

# a few functions to determine how many bits we in f.p. calculations
int bitTable[1024];
void initBitTable() {
  int i, j, k;
  bitTable[0] = 10;
  for (k = 0, j = 1; j < 1024; k++, j *=2)
    for (i = 0; i<1024; i += j)
      bitTable[i] = k;
}

const MM 0
const EE 1
typedef int meType[2];

# Compute # mantissa bits (omit trailing zeros) and the exponent
void meBits(double x, meType result) {
  int i;
  if (x < 0) x = -x;
  if (x == 0) {    result[EE] = 0; result[MM] = 0; }
  else {
    frexp(x, result+EE);
    i = bitTable[(int)fmod(x, 1024.0)];
    result[MM] = (i > result[EE])? 0 : result[EE]-i;
  }
}

const MAX(a,b) ((a)>(b) ? (a) : (b))

void meAdd(meType a, meType b, meType result) {
  int d = a[EE]-b[EE];
  if (a[MM] == 0) { # add 0 to b
    result[MM] = b[MM];
    result[EE] = b[EE];
  } else if (b[MM] == 0) { # add 0 to a
    result[MM] = a[MM];
    result[EE] = a[EE];
  } else if (d==0) {
    result[MM] = MAX(a[MM], b[MM])+1; # include carry if equal exponent
    result[EE] = a[EE]+1; # should really do this other times, too.
  }
  else if (d<0) {
    result[MM] = MAX(a[MM]-d, b[MM]); # align exponents & add/subtract
    result[EE] = b[EE];
  } else {
    result[MM] = MAX(b[MM]+d, a[MM]);
    result[EE] = a[EE];
  }
}

void meMult(meType a, meType b, meType result) {
  if ((a[MM] == 0) || (b[MM] == 0)) {
    result[MM] = 0;
    result[EE] = 0;
  }
  else {
    result[MM] = a[MM]+b[MM]-1;
    result[EE] = b[EE]+a[EE];
  }
}

double spdot(psphereType sp, ppointType pv, ppointType sv) {


  meType sx, sy, sz, sq, dx, dy, dz, dq;
  meType tx, ty, tz, tq, t, dd;

  meBits((sp)->x, sx);
  meBits((sp)->y, sy);
  meBits((sp)->z, sz);
  meBits((sp)->sq,sq);
  meBits((pv)->x-(sv)->x,  dx);
  meBits((pv)->y-(sv)->y,  dy);
  meBits((pv)->z-(sv)->z,  dz);
  meBits((pv)->sq-(sv)->sq,dq);

  meMult(sx, dx, tx);
  meMult(sy, dy, ty);
  meMult(sz, dz, tz);
  meMult(sq, dq, tq);

  meBits((sp)->x*((pv)->x-(sv)->x), t);
  if ((t[EE] != tx[EE]) ||(t[EE] != tx[EE]))
    printf("(%2d;%2d)=(%2d;%2d) = %2d;%2d x %2d;%2d\n",
	   t[MM], t[EE], tx[MM], tx[EE], sx[MM], sx[EE], dx[MM], dx[EE]);
  meBits((sp)->y*((pv)->y-(sv)->y), t);
  if ((t[EE] != ty[EE]) ||(t[EE] != ty[EE]))
    printf("(%2d;%2d)=(%2d;%2d) = %2d;%2d x %2d;%2d\n",
	   t[MM], t[EE], ty[MM], ty[EE], sy[MM], sy[EE], dy[MM], dy[EE]);
  meBits((sp)->z*((pv)->z-(sv)->z), t);
  if ((t[EE] != tz[EE]) ||(t[EE] != tz[EE]))
    printf("(%2d;%2d)=(%2d;%2d) = %2d;%2d x %2d;%2d\n",
	   t[MM], t[EE], tz[MM], tz[EE], sz[MM], sz[EE], dz[MM], dz[EE]);
  meBits((sp)->sq*((pv)->sq-(sv)->sq), t);
  if ((t[EE] != tq[EE]) ||(t[EE] != tq[EE]))
    printf("(%2d;%2d)=(%2d;%2d) = %2d;%2d x %2d;%2d\n",
	   t[MM], t[EE], tq[MM], tq[EE], sq[MM], sq[EE], dq[MM], dq[EE]);


  meAdd(tx, ty, t);
  meAdd(t, tz, t);
  meAdd(t, tq, t);

  double d = (sp)->x*((pv)->x-(sv)->x) +(sp)->y*((pv)->y-(sv)->y)
    +(sp)->z*((pv)->z-(sv)->z) +(sp)->sq*((pv)->sq-(sv)->sq);

  meBits(d, dd);

  if ( (t[MM] > 53) || ((t[MM] - dd[MM]) > 30) )
    {
    printf("(%2d;%2d)=(%2d;%2d) : (%2d;%2d %2d;%2d %2d;%2d %2d;%2d)\n                 *<%2d;%2d %2d;%2d %2d;%2d %2d;%2d>\n                 = %2d;%2d+%2d;%2d+%2d;%2d+%2d;%2d\n",
	   dd[MM], dd[EE], t[MM], t[EE],
	   sx[MM], sx[EE], sy[MM], sy[EE], sz[MM], sz[EE], sq[MM], sq[EE],
	   dx[MM], dx[EE], dy[MM], dy[EE], dz[MM], dz[EE], dq[MM], dq[EE],
	   tx[MM], tx[EE], ty[MM], ty[EE], tz[MM], tz[EE], tq[MM], tq[EE]);
    }
  return d;
}




## d3.h

#- d3 Delaunay triangulator in 3d             Jack Snoeyink Aug 2003
Implements Watson's incremental Delaunay, with sphere-based search and
simplical data structure storing vertices and opposite neighbors.
Handles degeneracies by perturbing points by increasing infinitesimals.
Guaranteed for integer coordinates of 10 bits. No flat simplices.
-#

#include <math.h>

#- d3.c Delaunay/Power diagram function
d3batch takes input vertices, and returns a (compact) corner table for Delaunay.
REQUIRES that the first point is at infinity,
and that the first 5 are in general position. (I should verify or relax this.)
Since all points have radii assigned already, it can produce power diagrams.
-#
void d3batch(ppointType vert, int nvert, # input vertices (pointType[]) & number of vertices
        pcornerType *result, int *ncorners); # output corner table

#- int d3initialize(const pd3stateType this, ppointType vertArray, int nvert)
initialize d3state this from vertArray.  Allocates memory for corners and spheres,
sets up the free list for tetrahedra, and creates spheres for the first five points.
REQUIRES: vertArray[0] contains the point at infinity and vertArray[1..4] contain
four finite points; these first five points must be in general position.
-#
int d3initialize(const pd3stateType this, ppointType vertArray, int nvert);

# LOCATION ROUTINES walk a mesh stored in s, active, & sph starting from corner start to find a point pv.
# result is the index of a simplex whose sphere strictly contains pv.
# (The simplex itself need not contain pv.)
# The return code is positive if we succeed (# of location steps at present)
# Return code of 0 means that we failed, perhaps because points with small weight have no Voronoi cells
# Return code of -1 means that we found the duplicate of a vertex in the mesh. IN THIS CASE, result
#   is the location of the corner where we found the duplicate!
int d3locSphere(const pd3stateType this, ppointType pv, int start,
                    int *result); # output

#- d3insert(const pd3stateType this, int vi, int p)
   inserts the point this->vert[vi] that is contained in sphere p
   into the delaunay triangulation stored in this.
   (d3locSphere may be used to obtain p.)
 -#
void d3insert(const pd3stateType this, int vi, int p);

#- d3compactCorners(const pd3stateType this, pcornerType *result, int *ncorners);
   Takes corner table this->s in which some corners/tetrahedra are unused, and
   returns a compactified corner table result, and its length ncorners.
   DESTROYS the corner table s and active flags in the process.
-#
int d3compactCorners(const pd3stateType this, # arrays this->s & this->active are DESTROYED!
                    pcornerType *result, int *ncorners); # output

# Functions to access corner tables
const MOD4(a) (a & 3)
const TETRA(corner) ((corner) >> 2)
const INDEX(corner) (MOD4(corner))
const CORNER(tetra,index) (((tetra)<<2)+(index))
const BASECORNER(corner) ((corner) & 0xFFFFFFFC)
const LASTCORNER(corner) ((corner)|3)

const DEAD(p) (this->active[p] <= 0) # is this a dead or killed tetrahedron?
const KILL(p) {this->active[p] = -1;} # kill tetrahedron
#const DEAD(p) (sph[p].sq < 0) # is this a dead tetrahedron?
#const KILL(p) {sph[p].sq = -1;} # kill tetrahedron

const infiniteV(pv,vert) ((pv) == vert) # first point is at infinity
const infiniteC(c) infiniteV(this->s[c].v, this->vert) # corner uses inf pt
const infiniteP(p) (this->sph[p].sq == 0.0) # if tetra uses infinite point

#- Tetrahedra are groups of four corners in order of increasing
vertex index, except that the first two may be swapped to ensure
that the orientation determinant is positive.
I.e., take the lex smallest alternating permutation with positive sign.
There must always be an odd number of swaps between two permutations;
we swap the first two if necessary to achieve this.
-#

# set corner's vertex and opposite in tetrahedron structure
#void setCornerVC(int c, int vv, int op) {
const setCornerVC(c, vv, op) { this->s[c].v = vv; this->s[c].opp = op;}
const setCornerPairV(c, op, vv, ov) { setCornerVC(c, vv, op); setCornerVC(op, ov, c); }
const setCornerVCN(c, vv, op) { setCornerVC(c, vv, op); this->s[op].opp = c; } # set corner & adjust nhbr opp

## d3.c

#- d3 Delaunay triangulation in 3d             Jack Snoeyink Aug 2003
   Implements Watson's incremental Delaunay, with sphere-based location,
   and simplical data structure storing vertices and opposite neighbors.
   Handles degeneracies by perturbing points by increasing infinitesimals.
   Guaranteed for integer coordinates of 10 bits. No flat simplices.
-#

#include "delaunay3.h"
#include "d3.h"


const spdot(sp,pv,sv) ((sp)->x*((pv)->x-(sv)->x)+(sp)->y*((pv)->y-(sv)->y)\
         +(sp)->z*((pv)->z-(sv)->z)+(sp)->sq*((pv)->sq-(sv)->sq))
const spdotInf(sp,pv) (sp)->sq


#- Allocate or free space for a tetrahedron.
   When allocating, this->liveTetra is the location fo the new tetrahedron.
-#
const STARTTETRA(this) \
{ startTetraCnt++; (this)->liveTetra = (this)->freeTetra; \
  (this)->freeTetra = (this)->s[CORNER((this)->liveTetra,0)].opp;  \
  ASSERT(DEAD((this)->liveTetra), "Reusing existing tetrahedron?"); (this)->active[(this)->liveTetra] = 1; \
  if ((this)->maxTetra <= (this)->liveTetra) { (this)->maxTetra = (this)->liveTetra+1; \
  if ((this)->maxTetra >= (this)->limitmaxTetra) { printf("AUDIT: %d > limitmaxTetra\n", (this)->liveTetra); exit(EXIT_FAILURE); }}#--#\
} # AUDIT

const FREETETRA(this, p) \
{ freeTetraCnt++;  #-TESTING-# ASSERT((this)->active[p]<0, "Freeing already free tetrahedron?"); #--#\
  (this)->active[p]=0;  if (((this)->maxTetra-p)<FORGETLIMIT) {\
  (this)->s[CORNER(p,0)].opp = (this)->freeTetra; (this)->freeTetra = p; #--#  \
}}

#Tetrahedron manipulation tables
const short offset[4][4] = { # c+offset[i][INDEX(c)] advances c to (c+i)mod5
  #-0-#{0,0,0,0}, #-1-#{1,1,1,-3}, #-2-#{2,2,-2,-2}, #-3-#{3,-1,-1,-1}};
const INCREMENT(c) (c+offset[1][INDEX(c)])
# drop[i] contains new vertex order after vertex i is dropped and replaced by pv on same side.
# offdr[i] contains drop(i)-index(i)
# invdrop[i][k] = j whenever drop[i][j] = k.  4s signal i=k; bad because i is dropped.
const short drop[4][3] = {{2,1,3}, {0,2,3}, {1,0,3}, {0,1,2}};
const short offdr[4][3] = {{2,1,3}, {-1,1,2}, {-1,-2,1}, {-3,-2,-1}};
const short invdrop[4][4] = {{4,1,0,2}, {0,4,1,2}, {1,0,4,2}, {0,1,2,4}};

#include "d3audit.c"

const double TWO25 = 1024.0*1024.0*32.0;

inline double fm(double x) {
  double y = fmod(x, TWO25);
  return y<0 ? y+TWO25: y;
}

inline double InSpherev(psphereType sp, ppointType pv, ppointType sv) {
  double d = spdot(sp, pv, sv); # Return true if inside (==negative)
#ifndef NOASSERT
  if ((fabs(d/((sp)->x)*((pv)->x-(sv)->x)) < 1e-10) # massive cancellation
      || (fabs(d/((sp)->y)*((pv)->y-(sv)->y)) < 1e-10)
      || (fabs(d/((sp)->z)*((pv)->z-(sv)->z)) < 1e-10)
      || (fabs(d/((sp)->sq)*((pv)->sq-(sv)->sq)) < 1e-10)) {

  double dmod = fm(fm((sp)->x)*((pv)->x-(sv)->x)+fm((sp)->y)*((pv)->y-(sv)->y)+fm((sp)->z)*((pv)->z-(sv)->z)+fm((sp)->sq)*((pv)->sq-(sv)->sq));
  #-  if (fm(d) != dmod) {
    printf("%% %10.0lf != %10.0lf: %18.0lf - %10.0lf = InSphere(<%8.0f %8.0f %8.0f %14.0f>*( %5.0f %5.0f %5.0f %10.0f))\n",
	   fm(d), dmod, d, fm(d)-dmod, sp->x,sp->y,sp->z,sp->sq, pv->x-(sv)->x, pv->y-(sv)->y,pv->z-(sv)->z,pv->sq-(sv)->sq);
  } #--#
    double d2 = d -(fm(d)-dmod);
    if ( (d<0 && d2 >=0) || (d==0 && d2 != 0) || (d>0 && d2 <=0)) {
      printf("%%%% ERROR: sgn(%lf) != sgn(%lf)\n", d, d2);
      #  if (fm(d) != dmod) {
    printf("%% %10.0lf != %10.0lf: %18.0lf - %10.0lf = InSphere(<%8.0f %8.0f %8.0f %14.0f>*( %5.0f %5.0f %5.0f %10.0f))\n",
	   fm(d), dmod, d, fm(d)-dmod, sp->x,sp->y,sp->z,sp->sq, pv->x-(sv)->x, pv->y-(sv)->y,pv->z-(sv)->z,pv->sq-(sv)->sq);
    } }
#endif
  inSphereCnt++; #-TESTING-#
  return d; # perturb those on sphere to inside

}

inline void makeSphereV(psphereType sp, ppointType v0, ppointType v1, ppointType v2, ppointType pv, ppointType vert) {
  double x0, y0, z0, sq0, x1, y1, z1, sq1, x2, y2, z2, sq2;
  double xy, xz, xs, yz, ys, zs; # 2x2 minors
  # make sphere: only v0 or v1 may be infinte.
  sphereCnt++; #-TESTING-#
  if(!infiniteV(v0, vert)) {
    x0 = v0->x - pv->x; y0 = v0->y - pv->y;
    z0 = v0->z - pv->z; sq0 = v0->sq - pv->sq;
  } else { x0 = v0->x; y0 = v0->y; z0 = v0->z; sq0 = v0->sq; }
  if(!infiniteV(v1, vert)) {
    x1 = v1->x - pv->x; y1 = v1->y - pv->y;
    z1 = v1->z - pv->z; sq1 = v1->sq - pv->sq;
  } else { x1 = v1->x; y1 = v1->y; z1 = v1->z; sq1 = v1->sq; }
  x2 = v2->x - pv->x; y2 = v2->y - pv->y;
  z2 = v2->z - pv->z; sq2 = v2->sq - pv->sq;
  xy = DET2(0,1,x,y);
  xz = DET2(0,1,x,z);
  yz = DET2(0,1,y,z);
  xs = DET2(0,1,x,sq);
  ys = DET2(0,1,y,sq);
  zs = DET2(0,1,z,sq);
  sp->x  = -y2*zs +z2*ys -sq2*yz;
  sp->y  =  x2*zs -z2*xs +sq2*xz;
  sp->z  = -x2*ys +y2*xs -sq2*xy;
  sp->sq =  x2*yz -y2*xz +z2*xy;
  #  sp->w  = -p2->x*sp->x -p2->y*sp->y -p2->z*sp->z -p2->sq*sp->sq;
  #-  printf("disp('Sphere equation: <%5.0f %5.0f %5.0f %5.0f %5.0f>')\n",sp->x,sp->y,sp->z,sp->sq);
      (void)fflush(stdout);
      printf("%%Sp ck: %g %g %g %g\n", spdot(sp,v0), spdot(sp,v1), spdot(sp,v2), spdot(sp,pv)); #--#
}


# when we initialize, this is what we fill in.
#const int initialopp[] = {11,5,15,21,25, 1,10,16,20,26, 6,0,17,22,27, 2,7,12,23,28, 8,3,13,18,29, 4,9,14,19,24};
const int initialopp[5][4] = {
  {CORNER(1,1), CORNER(2,0), CORNER(3,1), CORNER(4,0)},
  {CORNER(2,1), CORNER(0,0), CORNER(3,0), CORNER(4,1)},
  {CORNER(0,1), CORNER(1,0), CORNER(3,2), CORNER(4,2)},
  {CORNER(1,2), CORNER(0,2), CORNER(2,2), CORNER(4,3)},
  {CORNER(0,3), CORNER(1,3), CORNER(2,3), CORNER(3,3)}};

const int initialv[2][5][4] = {{{1,2,3,4}, {2,0,3,4}, {0,1,3,4}, {1,0,2,4}, {0,1,2,3}},
                                      {{0,2,3,4}, {2,1,3,4}, {1,0,3,4}, {0,1,2,4}, {1,0,2,3}}};


#- void d3initialize(const pd3stateType this, ppointType vertArray, int nvert)
initialize d3state this from vertArray.  Allocates memory for corners and spheres,
sets up the free list for tetrahedra, and creates spheres for the first five points.
REQUIRES: vertArray[0] contains the point at infinity and vertArray[1..4] contain
four finite points; these first five points must be in general position.
It permutes to get this to happen
-#
int d3initialize(const pd3stateType this, ppointType vertArray, int nvert)
{
  int itry, j, p, last;
  double d;
  pointType vtemp;

  #  initBitTable(); # BITS
  this->vert = vertArray;
  this->limitmaxTetra = TETperV*nvert; # allocate space for spheres and corners
  this->s   = (pcornerType) calloc(4*this->limitmaxTetra, sizeof(cornerType)); # per corner: v, opp
  this->sph = (psphereType) calloc(this->limitmaxTetra, sizeof(sphereType)); # per tetra: sphere eqn
  this->active = (int *) calloc(this->limitmaxTetra, sizeof(int)); # flag -1 unused, 0 dead, 1 alive

  if ((this->s == NULL) || (this->sph == NULL) || (this->active == NULL)) {
    printf("ERROR: d3batch could not calloc memory for data structures\n");
    exit(EXIT_FAILURE);
  }

  # initialize tetrahedra
  last = -1; #- set up free list of tetrahedra -#
  this->freeTetra = this->limitmaxTetra;
  do {
    this->freeTetra--;
    this->active[this->freeTetra] = 0; # KILL(freeTetra);
    this->s[CORNER(this->freeTetra,0)].opp = last;
    last = this->freeTetra;
  } while (this->freeTetra > 5);

  this->active[4] = 2; # create first sphere
  itry = 5; j = 0; # permute if first five are not in general position (Thanks, Leo)
  makeSphereV(this->sph+4, this->vert+0, this->vert+1, this->vert+2, this->vert+3, this->vert);
  d = spdot(this->sph+4, this->vert+4, this->vert+3); # if d<0, then we need to swap

  while (d == 0.0)
    {
      #	printf("Permute points to get first five into general position\n");
      if (itry>=nvert) # failure
	return 0;
      j = j%4+1;
      vtemp = *(vertArray+j);
      *(vertArray+j) = *(vertArray+itry);
      *(vertArray+itry) = vtemp;
      itry++;
      makeSphereV(this->sph+4, this->vert+0, this->vert+1, this->vert+2, this->vert+3, this->vert);
      d = spdot(this->sph+4, this->vert+4, this->vert+3); # if d<0, then we need to swap
    }

  if (d < 0) {
    this->sph[4].x = -this->sph[4].x; this->sph[4].y = -this->sph[4].y; this->sph[4].z = -this->sph[4].z;
    this->sph[4].sq = -this->sph[4].sq;
  }

  for (p=0; p<5; p++) {
    for (j=0; j<4; j++) {       # pay attention to orientation when assigning vertices
      setCornerVC(CORNER(p,j), this->vert+initialv[d<0][p][j], initialopp[p][j]);
    }
    makeSphereV(this->sph+p, this->s[CORNER(p,0)].v, this->s[CORNER(p,1)].v, this->s[CORNER(p,2)].v, this->s[CORNER(p,3)].v, this->vert);
    this->active[p] = 1;
    # swap first two if d<0.
    if (((d<0)&&(p==1)) || ((d>0) && (p==0)))
      {ASSERT(this->sph[p].sq >0, "Somehow vertp at infinity is in or on sphere p in init");}
    else
      {ASSERT(spdot(this->sph+p, this->vert+p + (d<0)*(p<2)*(1-2*p), this->s[CORNER(p,3)].v) > 0,
              "Somehow vertp is in or on sphere p in init.");}
  }
  this->liveTetra = 4;
  this->maxTetra = 5;

  return 1; # success
}



const AUDITDROP dropCnt++;#printf("Dropping pv %d = (%5d; %5.0f %5.0f %5.0f; %10f)\n", pv-this->vert, pv->index, pv->x, pv->y, pv->z, pv->sq);  cornerPrint4(this, *result)

# LOCATION ROUTINES walk a mesh stored in s, active, & sph starting from corner start to find a point pv.
# result is the index of a simplex whose sphere strictly contains pv.
# (The simplex itself need not contain pv.)
# The return code is positive if we succeed (# of location steps at present)
# Return code of 0 means that we failed, perhaps because points with small weight have no Voronoi cells
# Return code of -1 means that we found the duplicate of a vertex in the mesh. IN THIS CASE, result
#   is the location of the corner where we found the duplicate!
int d3locSphere(const pd3stateType this, ppointType pv, int start,
                    int *result) { # output
  int j, guard; # loop variables
  int c1, c2; # corners
  psphereType s1, s2; # spheres
  double I1, I2, d; # Insphere values

  #  printf("LOC:  %d = (%5d; %5.0f %5.0f %5.0f; %10f)\n", pv-this->vert, pv->index, pv->x, pv->y, pv->z, pv->sq);

  guard = 2*this->maxTetra+4; # prevent infinite loops
  c1 = CORNER(start,0); # corner in start
  s1 = this->sph+start; # sphere at start
  I1 = InSpherev(s1, pv, this->s[c1+3].v); # Check if strictly inside start sphere.
  if (I1 < 0) {# found already
    *result = start;
    return 1; # success on first try
  }

  while (--guard) {
    c2 = this->s[c1].opp; s2 = this->sph + TETRA(c2); # nhbr corner, sphere, value
    I2 = InSpherev(s2, pv, this->s[LASTCORNER(c2)].v);
    if (I2 < 0) {  # found one!
      *result = s2 - this->sph;
      return 2*this->maxTetra+5-guard; # number of steps
    }
    d = s2->sq * I1 - s1->sq * I2; # Warning: if s1 & s2 are same sphere, this is zero
    locateSideCnt++; # STATS
    if (d==0) { # We may be on two spheres---check for duplicate vertex
      if (EQUALPV(pv,this->s[c2].v)) { *result = c2; AUDITDROP; return -1; }
      j = INDEX(c2);
      if (EQUALPV(pv,this->s[c2+offset[1][j]].v)) { *result = c2+offset[1][j]; AUDITDROP; return -1; }
      if (EQUALPV(pv,this->s[c2+offset[2][j]].v)) { *result = c2+offset[2][j]; AUDITDROP; return -1; }
      if (EQUALPV(pv,this->s[c2+offset[3][j]].v)) { *result = c2+offset[3][j]; AUDITDROP; return -1; }
      # otherwise no duplicate; we probably have s1 == s2. (Rare in protein data.)
      d = RANDBIT-0.5; # choose a random direction, as a hack. (Should do plane computation)
    }
    if (d < 0) # if on I1 side
      c1 = INCREMENT(c1);
    else {# on I2 side
      c1 = INCREMENT(c2); s1 = s2; I1 = I2;
    }
  }
  result = 0; # location failure
  return 0;
}

#- d3compactCorners(const pd3stateType this, pcornerType *result, int *ncorners);
   Takes corner table this->s in which some corners/tetrahedra are unused, and
   returns a compactified corner table result, and its length ncorners.
   DESTROYS the corner table s and active flags in the process.
   returns false if it is unable to allocate memory.
-#
int d3compactCorners(const pd3stateType this, # arrays this->s & this->active DESTROYED!
                    pcornerType *result, int *ncorners) { # output
  int c, nc, i, j;
  pcornerType pc;

  i = 0; # count actives
  for (j = 0; j < this->maxTetra; j++)
    this->active[j] = (DEAD(j))? -1 : i++; # make old->new pointer dictionary for active tetra

  *ncorners = 4*i; # number of corners to return
  *result = pc = (pcornerType) calloc(*ncorners, sizeof(cornerType)); # return corner table: v, opp
  if (pc != NULL) {
    for (j = 0; j < this->maxTetra; j++) # compact the corners
      if (this->active[j] >= 0) {
        c = CORNER(j,0);  # old corner c -->  new corner ptr pc = (*result)+0,1,2,...
        pc->v = this->s[c].v; nc = this->s[c].opp;
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc)); # assign new tetra # w/ old index
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp;
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp;
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++; c++;
        pc->v = this->s[c].v; nc = this->s[c].opp;
        pc->opp = CORNER(this->active[TETRA(nc)], INDEX(nc));
        pc++;
      }
  }
  free(this->s);  # free old corner list
  free(this->active);
  return (pc != NULL);  # true if we were successful
}



#- d3insert(const pd3stateType this, int vi, int p)
   inserts the point this->vert[vi] that is contained in sphere p
   into the delaunay triangulation stored in this.
   (d3locSphere may be used to obtain p.)
 -#
void d3insert(const pd3stateType this, int vi, int p) {
  int i, j, off;
  int b, c, newb; # corners
  int nc, ni, dead, jdead; # indices
  double d;
  ppointType v0, v1, v2; # pointers to vertices
  ppointType pv = this->vert+vi;

    # Tetrahedra containing pv are "dead", and are pushed onto kill stack.
    # We use DFS with stack pst to find them and kill them
    # At live-dead boundary, we save dead tetras on stack nhbr,
    #  then make new tetras and hook in to live by setting the last opp pointer.
    #
    # Invariants/operations: Tetrahedron p is marked alive or dead on first visit.
    #    Corner c is pushed on stack when TETRA(this->s[c].opp) is marked dead.
    #
    # On termination, stack nhbr contains dead corners with live neighbors
    #    that have new tetras (so this->s[nhbr].opp != this->s[this->s[nhbr].opp].opp temporarily.)
    #    Stack kill contains old tetrahedra for final recycling.

    stkINIT(this->dfs);  # DFS stack holds corners opposite dead tetras
    stkINIT(this->idfs); # iDFS stack holds corners opposite infinite tetras with pv on bdry
                         #    (these are a special case: dead, but don't propagate)
    stkINIT(this->nhbr); # stack for dead corners with live nhbr tetras
    stkINIT(this->kill); # stack of dead tetras to recycle
    b = CORNER(p,0);
    PUSH(p, this->kill); KILL(p); # kill tetra initial p,
    PUSH(this->s[b++].opp, this->dfs); # stack neighbors
    PUSH(this->s[b++].opp, this->dfs);
    PUSH(this->s[b++].opp, this->dfs);
    PUSH(this->s[b  ].opp, this->dfs);

    while (!isEMPTY(this->dfs)) {
      c = POP(this->dfs); p = TETRA(c);
      #-        printf("::Popping %d with opp %d \n", c, this->s[c].opp);#--#
      ASSERT(DEAD(TETRA(this->s[c].opp)), "dfs stack element with non-dead neighbor");
      if (DEAD(p)) continue; # dead already
      d = InSpherev(this->sph + p, pv, this->s[LASTCORNER(c)].v); # Is pv in, out, or on?
      if (d < 0) { # kill and continue dfs if pv is strictly inside
        KILL(p); PUSH(p, this->kill); # kill and stack tetra
        j = INDEX(c);
        PUSH(this->s[c+offset[1][j]].opp, this->dfs); # stack neighbors to check
        PUSH(this->s[c+offset[2][j]].opp, this->dfs);
        PUSH(this->s[c+offset[3][j]].opp, this->dfs);
      }
      else if (d > 0  || this->sph[p].sq > 0) {  # pv is outside (or on with sp finite), so
	                                         # c is live neighbor of dead opp tetra this->s[c].opp
        PUSH(this->s[c].opp, this->nhbr); # remember old corner, so we can hook tetra into mesh later
        STARTTETRA(this); # make new tetrahedron liveTetra
        newb = CORNER(this->liveTetra,3); # last corner of new tetra
        setCornerVCN(newb, pv, c); # last corner is pv; also set opposite corner c. Do rest later.
      }
      else { # d==0 && sph[p] is infinite: handle two special cases
          if (this->sph[TETRA(this->s[c].opp)].sq == 0) { # if dead sphere is infinite, too
            PUSH(c, this->idfs); # then if c stays alive, we make tetra to it (flat, but infinite).
          } else { # dead sphere is finite; kill c and make tetras to neighbors, if they stay alive.
            KILL(p); PUSH(p, this->kill); # kill and stack tetra
            j = INDEX(c);
            PUSH(this->s[c+offset[1][j]].opp, this->idfs); # stack neighbors to check
            PUSH(this->s[c+offset[2][j]].opp, this->idfs);
            PUSH(this->s[c+offset[3][j]].opp, this->idfs);
          }
        }
    }

    while (!isEMPTY(this->idfs)) { # check the neighbors of infinite tetrahedra
      c = POP(this->idfs); p = TETRA(c);
      #-        printf("::Popping %d with opp %d \n", c, this->s[c].opp);#--#
      ASSERT(DEAD(TETRA(this->s[c].opp)), "dfs stack element with non-dead neighbor");
      if (DEAD(p)) continue; # dead already
      ASSERT(DEAD(TETRA(this->s[c].opp)), "Live corner c should have dead neighbor");
      PUSH(this->s[c].opp, this->nhbr); # remember old corner, so we can hook tetra into mesh later
      STARTTETRA(this); # make new tetrahedron liveTetra
      newb = CORNER(this->liveTetra,3); # last corner of new tetra
      setCornerVCN(newb, pv, c); # last corner is pv; also set opposite corner c. Do rest later.
    }

    # Now, we have stack of dead neighbors of live tetras, and we've hooked new tetras to them.
    while (!isEMPTY(this->nhbr)) {
      dead = POP(this->nhbr); jdead = INDEX(dead); #  dead tetra and index of dropped corner.
      #-        printf("--Popped %d(%d)\n", dead, jdead); #--#
      ASSERT(DEAD(TETRA(dead)), "corner on nhbr stack is not dead!?");
      newb = this->s[this->s[dead].opp].opp-3; # base of new tetra.

      dead -= jdead; # just use base of dead one.
      # new tetra has 0,1,2,3=pv;
      # corresponding old indices before jdead is dropped:
      #   drop[j][0],..,drop[j][3], (no corresp to pv)
      j = jdead;
      i = drop[jdead][0]; # old index of new corner 0;
      c = dead+i; # note i = INDEX(C);
      v0 = this->s[c].v; # copy vertex v0
      nc = this->s[c].opp; # go to neighbor
      # In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      # To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) {
        ni = INDEX(nc); off = indoff[i][ni][j]; # where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; # fix new j, i, c, and try neighbor
      }
      nc = this->s[nc].opp; # go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v0, nc-3+invdrop[i][j]); newb++;

      j = jdead;
      i = drop[jdead][1]; # old index of new corner 1;
      c = dead+i; # note i = INDEX(C);
      v1 = this->s[c].v; # copy vertex v1
      nc = this->s[c].opp; # go to neighbor
      # In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      # To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) {
        ni = INDEX(nc); off = indoff[i][ni][j]; # where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; # fix new j, i, c, and try neighbor
      }
      nc = this->s[nc].opp; # go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v1, nc-3+invdrop[i][j]); newb++;

      j = jdead;
      i = drop[jdead][2]; # old index of new corner 2;
      c = dead+i; # note i = INDEX(C);
      v2 = this->s[c].v; # copy vertex v2
      nc = this->s[c].opp; # go to neighbor
      # In tetra opp c, find new location of j.  That is new c. New j = INDEX(c.opp).
      # To avoid index calculations, maintain i = INDEX(c), nc = this->s[c].opp, ni = INDEX(nc).
      while (DEAD(TETRA(nc))) {
        ni = INDEX(nc); off = indoff[i][ni][j]; # where j goes relative to i is our new i.
        j = ni; i = ni + off; c = nc + off; nc = this->s[c].opp; # fix new j, i, c, and try neighbor
      }
      nc = this->s[nc].opp; # go to new tetra
      ASSERT(this->s[nc].v == pv, "Expected to find new tetra using pv after walking dead tetras. ");

      setCornerVC(newb, v2, nc-3+invdrop[i][j]); newb++;

      c = this->s[this->s[dead+jdead].opp].opp;
      ASSERT(v0==this->s[CORNERINOPP(0,c)].v, "v0 does not line up");
      ASSERT(v1==this->s[CORNERINOPP(1,c)].v, "v1 does not line up");
      ASSERT(v2==this->s[CORNERINOPP(2,c)].v, "v2 does not line up");

      makeSphereV(this->sph+TETRA(newb), v0, v1, v2, pv, this->vert); # use either this or makeSphereP above
    }
}



#- d3.c Delaunay/Power diagram function
   d3batch takes lifted input vertices, and returns a (compact) corner table for Delaunay.
   REQUIRES that the first point is at infinity,
   and that the first 5 are in general position. (I should verify or relax this.)
   Since all points have radii assigned already, it can compute power diagrams.
-#
void d3batch(ppointType vertArray, int nvert, # input vertices (pointType[]) & number of vertices
        pcornerType *result, int *ncorners)  # output corner table
{
  int k;
  int p; # tetra
  int vi; # vertex index in outermost loop
  ppointType pv; # vertex pointer in outermost loop
  pcornerType pc; # corner for copying to result
  d3stateType d3; # state data structure
  const pd3stateType this = &d3;

  if (!d3initialize(this, vertArray, nvert)) {
    printf("ERROR: Can't initialize\n");
    exit(EXIT_FAILURE);
  }
#ifdef STATS
  this->killmax = this->dfsmax = this->idfsmax = this->nhbrmax = -1; # init STATS
#endif
  maxLocate = locateSideCnt = sphereCnt = startTetraCnt = freeTetraCnt= inSphereCnt=0;

  for (vi = 5; vi < nvert; vi++) { # incrementally insert vert[vi]
    #LOCATE: find some tetrahedron with sphere strictly containing vert[vi]
    if ((k = d3locSphere(this, this->vert+vi, this->liveTetra, &p)) < 1)
      continue; # if we fail to locate (duplication or other reason) just skip vert[vi]
    if (k > maxLocate) maxLocate = k; # STATS
#ifndef NOASSERT
    if (vi%2000 == 0) { # print a few search results
      printf("%%Found %3d in %4d(%3d) after %d steps %f\n",
             vi, CORNER(p,0), p, k, InSpherev(this->sph+p, this->vert+vi, this->s[CORNER(p,3)].v));
      auditCorners(this, 1); # audit, and check spheres, too
      (void)fflush(stdout);
    } #--#
#endif
    d3insert(this, vi, p); #  insert vertex vi, which is in sphere of tetra p.
    while (!isEMPTY(this->kill)) { # recycle memory of dead tetrahedra
      p = POP(this->kill);
      FREETETRA(this, p);
    }
  }
#ifdef STATS
  printf(" %10d\tPoints given\n %10d\tPoints used\n %10d\t Max Tetra\n", vi, vi-dropCnt, this->maxTetra);
      printf("We performed:\n %10d\tinSphere tests\n %10d\tplane tests\n %10d\ttetrahedra created-\n %10d\tfreed =\n %10d\ttetrahedra\n %10d\tsphere equations computed.\n",
             inSphereCnt, locateSideCnt, startTetraCnt, freeTetraCnt, startTetraCnt-freeTetraCnt, sphereCnt);
      printf("Max locate steps=%d, max killed=%d, max inf=%d, max created=%d, maxdfs=%d\n",
	     maxLocate, this->killmax, this->idfsmax, this->nhbrmax, this->dfsmax);
      fflush(stdout);
      auditCornersAux(this,1); #-TESTING-#
#endif
  auditCornersAux(this,1); #-TESTING-#

  free(this->sph);  # done with spheres; free them

  if (!d3compactCorners(this, result, ncorners)) { # disposes s and active!!
    printf("ERROR: d3.c could not allocate corner table to return\n");
    exit(EXIT_FAILURE);
  }
}

## delaunay3.h

#- delaunay3 Delaunay/Power diagrams in 3d Jack Snoeyink Aug 2003
Implements Watson's incremental Delaunay, with sphere-based search and
simplical data structure storing vertices and opposite neighbors.
Points must be scaled to integers; guaranteed for differences of 10
bits, and uses leveling and hilbert curve to guarantee 16 bits if the
points are well distributed.  Works for pdb files, which are 20 bits.
Handles degeneracies by perturbing points by increasing infinitesimals
to guarantee that all simplices are full-dimensional.
-#

#ifndef DELAUNAY3_H
const DELAUNAY3_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef float coord; # We'll just use floats for x,y,z coordinates

typedef struct { # Basic point data structure element
  int index; # index back to original line in file.
  coord x,y,z; double sq;  # coordinates: 3 spatial + lifted
} pointType, *ppointType;

#tetrahedra are groups of four consecutive corners
typedef struct cornerType { # data associated with each corner
  ppointType v; # index of vertex
  int opp; # pointer to opposite corner in neighboring tetra
} cornerType, *pcornerType;

#- Tetrahedra are groups of four corners in order of increasing
vertex index, except that the first two may be swapped to ensure
that the orientation determinant is positive.
I.e., take the lex smallest alternating permutation with positive sign.
There must always be an odd number of swaps between two permutations;
we swap the first two if necessary to achieve this.  -#

#- d3.c Delaunay/Power diagram function
d3batch takes input vertices, and returns a (compact) corner table for Delaunay.
REQUIRES that the first point is at infinity,
and that the first 5 are in general position.
(I should probably verify or relax this.)
Since all points have radii assigned already, it can handle power diagrams.
-#
void d3batch(ppointType vert, int nvert, # input vertices (pointType[]) & number of vertices
        pcornerType *result, int *ncorners); # output corner table


const XX 0
const YY 1
const ZZ 2

const NVERT     5000 # number of random vertices
const MAXVERT 400000 # pdb has 5 digit atom# field
const TETperV    10 #
const FORGETLIMIT  10000 # forget tetra after this many (for better locality of ref? Doesn't help.)

# Some useful definitions
const EQUALV(i,j) (vert[i].x == vert[j].x && vert[i].y == vert[j].y && vert[i].z == vert[j].z)
const EQUALPV(pv,v) (pv->x == v->x && pv->y == v->y && pv->z == v->z)

const DET2(p,q,i,j)((i##p)*(j##q) - (j##p)*(i##q))
const DOT(a,b) ((a)->x*(b)->x + (a)->y*(b)->y + (a)->z*(b)->z)

# Stack data structure operations
const STACKMAX 2000
const POP(stack) (stack##st[stack##sp--])
const isEMPTY(stack) (stack##sp < 0)
const stkINIT(stack) {stack##sp = -1; }

#ifndef STATS
const PUSH(value, stack) { stack##st[++stack##sp] = value; }
const stkDECLARE(stack,stn) int stack##sp, stack##st[STACKMAX];
#else
const PUSH(value, stack) { \
    stack##st[++stack##sp] = value; \
    if (stack##max < stack##sp) { stack##max = stack##sp; \
      if (stack##max >= STACKMAX) { \
        printf("ERROR: overflow stack %x pushing %d",  stack##st, value); exit(EXIT_FAILURE); } } #--#\
    }
const stkDECLARE(stack,stn) int stack##sp, stack##st[STACKMAX]; int stack##max; #AUDIT #--#
#endif

typedef struct sphereType { # sphere equation
  double x, y, z, sq; # Invariant sq > 0 for all created tetra, unless they use pt at infty
 } sphereType, *psphereType;

typedef struct {
  ppointType vert; # vertices: 0th is point at infinity!!!!
  pcornerType s;  # corner table
  psphereType sph; # spheres
  int *active; # which spheres are active
  int freeTetra; # head for free list for tetrahedra kept in opp[CORNER(tetra,0)]
  int liveTetra; # latest tetra; known to be live.
  int maxTetra; # AUDIT only
  int limitmaxTetra; # limit on # of created tetrahedra, spheres & corners/4.
  # stacks used in inserting pv
  stkDECLARE(dfs, "dfs");  # DFS stack to find dead tetras
  stkDECLARE(idfs, "idfs"); # DFS stack for tetras adj to infinite vertex (>4*30)
  stkDECLARE(nhbr, "nhbr"); # stack for dead corners with live neighbors
  stkDECLARE(kill, "kill"); # stack for base corners of tetras to recycle
} d3stateType, *pd3stateType;

# readers call setVert when they make a point
void setVert(ppointType v, int indx, double xx, double yy, double zz,
             double rad, double mult);

const VSUBASSN(vin, xx,yy,zz, vout) {\
    vout->index = vin->index;  vout->sq = vin->sq;\
    vout->x = vin->x - xx; vout->y = vin->y - yy; vout->z = vin->z - zz;\
  }

#Some compiler/unix variants
#ifndef FALSE
const FALSE 0
#endif

#ifndef TRUE
const TRUE 1
#endif

#ifndef strcasecmp
const strcasecmp(s1,s2)      strcmp(s1,s2)
const strncasecmp(s1,s2,n)   strncmp(s1,s2,n)
#endif

const HIGHCOORD 0x4000
const COORDMASK 0x3fff

#ifdef bcc
const RANDBIT (random(2))             # one random bit
const RAND2BIT (random(4))            # two random bits
const RANDPROB(mask) (random(mask+1)) # random bits masked
const RANDCOORD (double)random(HIGHCOORD);
const RANDOM(k) random(k)             # random number 0..k-1
#else
const RANDBIT (random()&1)            # one random bit
const RAND2BIT (random()&3)           # two random bits
const RANDPROB(mask) (mask&random())  # random bits masked
const RANDCOORD (double)(random()&COORDMASK)
const RANDOM(k) random()%(k)          # random number 0..k-1
#endif

#  fprintf(stderr, "\nASSERT FAILED (line %d of %s ): %s\n",
#ifndef NOASSERT
const ASSERT(bool, string) if (!(bool)) {\
  printf("\nASSERT FAILED (line %d of %s ): %s\n", \
          __LINE__, __FILE__, string); }
#else
const ASSERT(bool, string)
#endif

#endif

##delaunay3.c
#- delaunay3.c                   Jack Snoeyink Aug 2003
Read in a file (often pdb) and save index & coordinates
in a hierarchy file that is ordered by number of trailing zeros,
and a Hilbert curve on high bits.
This format is ready for Delaunay triangulation/power diagram
computation: bbox, num,
Flags:
-a asymmetric Hilbert curve
-b# # of bits each coord contribs to boxno
-l double levels
-g Gray code order
-z Z (Morton) order
-#

bboxType bb; # bounding box

int curvetype = 0, lflag = 0, nboxbits = 6;

#include "delaunay3.h"
#include "plyreader.h"
#include "pdbreader.h"

void BBprint() {
  printf("BB %d (%f %f %f; %f %f %f)\n", BBbits(bb),
	 bb[0], bb[1], bb[2], bb[3], bb[4], bb[5]);
}

const int boxshiftLUT[] = {0,3,6,9,12,15,18,21,24}; # shift boxno to get memory loc. 5->15

const int high[64] = {0, 1, 2,2, 3,3,3,3, 4,4,4,4,4,4,4,4,
		      5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};  # for BBbits

inline int BBbits(bboxType bb) { # returns the number of bits in bb range
  int range;
  int mask = 0xffffFFE0; # mask last 5 bits
  int shift = 0;
  range = bb[1]-bb[0]; # find max range in all three coords
  if (range < bb[3]-bb[2]) range = bb[3]-bb[2];
  if (range < bb[5]-bb[4]) range = bb[5]-bb[4];

  while (range & mask) { # some ones bits left above last 5.
    mask <<= 5;
    shift += 5;
  }
  return shift+high[range>>shift]; # shift needed to make all bits zero
}

#- for qsort: compares box numbers -#
int boxnoCmp(const void *p, const void *q) {
  return ((pointType *)p)->boxno - ((pointType *)q)->boxno;
}

#- addVB(b, vert) adds the vertex vert to the vertex buffer b.
Allocates memory for b->v if vert is the first, and
writes b->v to b->fname if b->v is full.
Increments total b->count, and  # in buffer, b->bcount.
-#
void addVB(pvbufferType b, ppointType vert) {
  int k;

  if (b->v == NULL)
    b->v = (ppointType) calloc(MAXBFSIZE, sizeof(pointType));
  else if (b->bcount == MAXBFSIZE) {
    if (b->fid == NULL) {
      b->fid = fopen(b->fname, "w");
      if (b->fid == NULL)
	printf("Null file in fopen %s\n", b->fname);
    }
    k = fwrite((void *) (b->v), sizeof(pointType), b->bcount, b->fid); # write out points
    if (b->v->boxno >= 0) printf("  W%s,%d,%o", b->fname, k, b->v->boxno);
    if (k != b->bcount)
      printf("ERROR on write points  j%d, err%d, %d\n", b->bcount, k, ferror(b->fid));
    b->bcount = 0;
  }
  b->v[b->bcount++] = *vert;
  b->count++;
}

#- writeUnique(nvert, vert[], boxshift, outfile)
sorts the nvert points in vert[] by vert[].boxno, then
writes a maximal set of distinct points to outfile.
Assumes that points with same coordinates have same boxno.
Side effect: shifts boxno by boxshift, and frees vert[].
Returns the number of points written.
-#
int writeUnique(int nvert, pointType vert[], int boxshift, FILE *outfile) {
  int nv = nvert;
  int j, k;

  if (nvert == 0)
    return 0;

  printf("writeUnique %d %x %d\n", nvert, vert, boxshift); fflush(stdout); #--#
  qsort(vert, nv, sizeof(pointType), boxnoCmp);
  for (j = nv; j > 0; j--) { # scan for duplicates
    k = j+1;
    while ((k < nvert)
	   && ((vert[k].boxno == -1) || (vert[k].boxno == vert[j].boxno))) {
      if ((vert[k].boxno != -1) && EQUALV(j,k)) {
	#	    printf("Found duplicate %d %d\n", j,k);
	nv--;
	vert[k].boxno = -1;
      }
      k++;
    }
  }

  #  printf(":writeUnique %d, dups %d \n", nv, nvert-nv); fflush(stdout); #--#
  for (j = 0; j < nvert; j++) { # for boxno that are still >= 0, shift and write
    #-    printf("    %d(%o,%o); (%d %d %d)\n", j,
	  vert[j].boxno, vert[j].boxno>>boxshift, vert[j].x, vert[j].y, vert[j].z); #--#
    if (vert[j].boxno >= 0) {
      vert[j].boxno >>= boxshift;
      #-      fprintf(outfile, "%d %d %d %d %d\n", vert[j].boxno, vert[j].index,
	      vert[j].x, vert[j].y, vert[j].z); #--#
      k = fwrite((void *) (vert+j), sizeof(pointType), 1, outfile); # write out point
      if (k != 1) {
	printf("ERROR on write points  j%d, err%d, %d\n", j, k, ferror(outfile));
	nv--;
      }#--#
    }
  }
  free(vert);
  printf("::writeUnique %d, dups %d \n", nv, nvert-nv); fflush(stdout); #--#
  return nv;
}

#- writeOrder(b, boxshift, outfile)
Sorts points in vertex buffer b by boxno, and writes them to outfile.
Depending on the number of points in b,
it either writes them from memory (reading them first if necessary),
or splits into NUMBIGBOX=64 buffers based on MASK6 bits of boxno>>boxshift,
and recursively writes the order for these buffers.
Returns the number of points written.
Side effects: empties file b->fname and frees b->v.
TODO: figure out the right way to delete a file for the DELETEFILE macro in .h
-#
int writeOrder(pvbufferType b, int boxshift, FILE *outfile) {
  int j,k;

  if (b->count <= 0)
    return 0;
  printf("writeOrder %s #%d, shift %d, so ", b->fname, b->count, boxshift);
  fflush(stdout); #--#
  if (b->count <= MAXBFSIZE) { # all points are in b->v; just write
    return writeUnique(b->bcount, b->v, boxshift, outfile);
  }
  # read all points into b->v, then write
  if ((b->count < RADIXSORTLIMIT) || (boxshift < BOXBITLIMIT)) {
    fclose(b->fid);
    b->v = (ppointType) realloc(b->v, sizeof(pointType)*b->count);
    b->fid = fopen(b->fname, "r");
    if (b->fid == NULL) {
      printf("Could not open buffer file for read %s.\n", b->fname);
      exit(1);
    }
    b->bcount += fread(b->v+b->bcount, sizeof(pointType), b->count, b->fid);
    if (b->bcount != b->count)
      printf("Counts %d != %d from %s\n", b->count, b->bcount, b->fname);
    fclose(b->fid);
    DELETEFILE(b->fname);
    return writeUnique(b->bcount, b->v, boxshift, outfile);
  }
  else { # split based on boxno>>boxshift and writeOrder recursively
    vbufferType vb[NUMBIGBOX];
    char fname[] = "tmpxxxx~";

    printf("splitting %d\n", boxshift);
    fclose(b->fid);
    b->fid = fopen(b->fname, "r");
    if (b->fid == NULL) {
      printf("Could not reopen buffer file for read %s.\n", b->fname);
      exit(1);
    }
    for (j = 0; j < NUMBIGBOX; j++) { # make buffers
      fname[4] = 'A'+boxshift/3;
      fname[5] = 'A'+(j>>3);
      fname[6] = '0'+(j&7);
      initVB(vb+j, fname);
    }
    j = b->bcount;
    for (k = 0; k < b->count; k++) {
      if (j <= 0) {
	if (MAXBFSIZE != (j =
			  fread(b->v, sizeof(pointType), MAXBFSIZE, b->fid)))
	  printf("At %d of %d, got only %d from buffer file %s.\n",
		 k, b->count, j, b->fname);
      }
      #-           printf("addVB(%s %o) ind %d shift %d box %o (%d %d %d)\n",
	      vb[((b->v[j].boxno>>boxshift)&MASK6)].fname,
	      ((b->v[j].boxno>>boxshift)), b->v[j].index, boxshift,
	      b->v[j].boxno,b->v[j].x,b->v[j].y,b->v[j].z);
	      fflush(stdout);#--#
      addVB(vb+((b->v[j].boxno>>boxshift)&MASK6), b->v+j);
      j--;
    }
    fclose(b->fid);
    DELETEFILE(b->fname);

    k = 0; # write out (and count) vbuffers
    for (j = 0; j < NUMBIGBOX; j++)
      k += writeOrder(vb+j, boxshift-6, outfile);
    printf("Done splitting %d, shift %d\n", k, boxshift);
    return k;
  }
}

#usage message
void usage()
{
    printf("makehier infile.pdb outfile.hier (flags) converts files to binary hierarchy for Delaunay\n");
    printf("makehier -n outfile.hier generates n random points\n");
    printf("     file types: ascii .txt (default), ascii .pdb, ascii .ply, hier (.hier)\n");
    printf("     flags: -a asymmetric Hilbert curve order\n");
    printf("            -b# number of leading coordinate bits to determine box number\n");
    printf("            -g Gray code order\n");
    printf("            -het (.pdb) include HETATM records\n");
    printf("            -l double the number of levels\n");
    printf("            -m# multiplier used in scaling coordinates to integers\n");
    printf("            -z Morton or Z order\n");
    exit(1);
}


#- reads or creates data, assigns boxno according to a space-filling curve,
and writes out in binary file in hierarchy with each level ordered by boxno.
-#
int main(int argc, char** argv) {
  FILE *fid; # input file id
  FILE *outfile; # output file id
  ppointType p; # input points
  coordInfoType coordinfo; # coordinate information

  vbufferType vb[NLEVELS]; # buffers for levels
  char line[132];
  int i, j, k, n;
  int highbit, boxlimit, boxshift;
  int totalvert;

  FILE * (*opener)(char *, int, char **); # declare function pointers
  void (*reader)(FILE *, ppointType *, int *);


  opener = txtopener; # default format is text
  reader = txtreader;

  # parse general parameters
  for (i = 1; i < argc; i++) {
    if ((0 == strncasecmp(argv[i],"-a", 2)))
      curvetype = 1;
    if ((0 == strncasecmp(argv[i],"-g", 2)))
      curvetype = 2;
    if ((0 == strncasecmp(argv[i],"-z", 2)))
      curvetype = 3;
    else if ((0 == strncasecmp(argv[i],"-l", 2)))
      lflag = TRUE;
    else if ((0 == strncasecmp(argv[i],"-b", 2))) {
      if (1 != sscanf(argv[i], "-b%d", &nboxbits))
	printf("couldn't sscanf bit number -b#\n");
    }
    else if ((0 == strncasecmp(argv[i],"-rnd", 4))) {
      opener = rndopener;
      reader = rndreader;
      if (1 != sscanf(argv[i], "-rnd%d", &nrandvert)) #
	if ((i+1 < argc) && (1 != sscanf(argv[i+1], "%d", &nrandvert)))
	  i++; # consume separate argument for # of random points
	else
	  nrandvert = 5000; #default # of random points
    }
  if (strstr(argv[i],".ply")) { # read .ply file
    opener = plyopener;
    reader = plyreader;
  }
  if (strstr(argv[i],".pdb") || strstr(argv[1],".ent")) { # read .pdb file
    opener = pdbopener;
    reader = pdbreader;
  }
  if (strstr(argv[i],".hier")) { # read .hier file
    opener = hieropener;
    reader = hierreader;
  }

  fid = opener(argv[i], argc, argv); # open the file
  if (fid == NULL)
    printf("Couldn't open %s\n", argv[1]);
  else
    d3read(fid, reader, &v, &nvert, &coordinfo); # read vertices & return coordinfo

  highbit = BBbits(bb)-1; # Choose which bits to use to define box numbers
  boxlimit = highbit+1 - nboxbits;
  if (boxlimit < BOXBITLIMIT)
    boxlimit = (highbit>BOXBITLIMIT) ? BOXBITLIMIT : highbit;
  printf("high bit %d; last box bit %d\n", highbit, boxlimit);
  boxshift = boxshiftLUT[highbit-boxlimit];

  printf("Boxshift %d\n", boxshift);

  # open output file
  outfile = fopen(argv[2], "w");
  printf("Open out %s\n", argv[2]);#--#
  if (outfile == NULL) {
    printf("Could not open final output file %s.\n", argv[2]);
    exit(1);
  }

  fprintf(outfile, "%%%%"); # save command line as first line of file
  for (i = 0; i < argc; i++)
    fprintf(outfile, "%s ", argv[i]);
  fprintf(outfile, "\n");

  # read from tmp files
  totalvert = 0;
  for (i = NLEVELS-1; i >= 0; i--)
    if (vb[i].count > 0) {
      printf("Level %d has %d points\n", i, vb[i].count);
      fflush(stdout);
      { # writeOrder for each level, but shift bbox and assign boxno, as well.
	int j,k;
	pvbufferType b = vb+i; # shorthand.

	printf("wo %s, %d, shift %d, so ", b->fname, b->count, boxshift);
	fflush(stdout); #--#
	if (b->count <= MAXBFSIZE) { # all points in buffer, just write
	  for (j = 0; j < b->count; j++) {
	    adjustCoords(b->v+j, -(bb[3]+bb[0])/2, -(bb[4]+bb[1])/2, -(bb[5]+bb[2])/2); # coords from 0 to bbrange
	    setBoxno(b->v+j, highbit, boxlimit);
	    #-	    printf("setboxno%d ind %d shift %d box %o (%d %d %d)\n", j, b->v[j].index,
		   boxshift, b->v[j].boxno,b->v[j].x,b->v[j].y,b->v[j].z);
	    fflush(stdout);#--#
	  }
	  totalvert += writeUnique(b->bcount, b->v, boxshift, outfile);
	}
	# read all points into buffer, then write
	else if (b->count < RADIXSORTLIMIT) {
	  fclose(b->fid);
	  b->v = (ppointType) realloc(b->v, sizeof(pointType)*b->count);
	  b->fid = fopen(b->fname, "r");
	  if (b->fid == NULL) {
	    printf("Could not open buffer file for read %s.\n", b->fname);
	    exit(1);
	  }
	  b->bcount += fread(b->v+b->bcount, sizeof(pointType), b->count, b->fid);
	  if (b->bcount != b->count)
	    printf("Counts %d != %d from %s\n", b->count, b->bcount, b->fname);
	  fclose(b->fid);
	  DELETEFILE(b->fname);
	  for (j = 0; j < b->count; j++) {
	    adjustCoords(b->v+j, -bb[0], -bb[2], -bb[4]); # coords from 0 to bbrange
	    setBoxno(b->v+j, highbit, boxlimit);
	    #-	    printf("setboxno%d ind %d shift %d box %o (%d %d %d)\n", j, b->v[j].index,
		   boxshift, b->v[j].boxno,b->v[j].x,b->v[j].y,b->v[j].z);
	    fflush(stdout);#--#
	  }
	  totalvert += writeUnique(b->bcount, b->v, boxshift, outfile);
	}
	else { # split into buffers and writeOrder for each
	  vbufferType vb2[NUMBIGBOX];
	  char fname[] = "tmpxxx~";

	  printf("splitting2: shift %d \n", boxshift);
	  fclose(b->fid);
	  b->fid = fopen(b->fname, "r");
	  if (b->fid == NULL) {
	    printf("Could not reopen buffer file for read %s.\n", b->fname);
	    exit(1);
	  }
	  for (j = 0; j < NUMBIGBOX; j++) { # make buffers
	    fname[3] = 'A'+boxshift/3;
	    fname[4] = 'A'+(j>>3);
	    fname[5] = '0'+(j&7);
	    initVB(vb2+j, fname);
	  }
	  j = b->bcount;
	  for (k = 0; k < b->count; k++) {
	    if (j <= 0) {
	      if (MAXBFSIZE != (j =
				fread(b->v, sizeof(pointType), MAXBFSIZE, b->fid)))
		printf("At %d of %d, got only %d from buffer file %s.\n",
		       k, b->count, j, b->fname);
	    }
	    adjustCoords(b->v+j, -bb[0], -bb[2], -bb[4]); # coords from 0 to bbrange
	    setBoxno(b->v+j, highbit, boxlimit);
	    #-	          printf("addVB(%s %o) ind %d shift %d box %o (%d %d %d)\n",
		    vb2[((b->v[j].boxno>>boxshift)&MASK6)].fname,
		    ((b->v[j].boxno>>boxshift)), b->v[j].index, boxshift,
		    b->v[j].boxno,b->v[j].x,b->v[j].y,b->v[j].z);
		    fflush(stdout);#--#
	    addVB(vb2+((b->v[j].boxno>>boxshift)&MASK6), b->v+j);
	    j--;
	  }
	  fclose(b->fid);
	  DELETEFILE(b->fname);

	  for (j = 0; j < NUMBIGBOX; j++)
	    totalvert += writeOrder(vb2+j, boxshift-6, outfile);
	}
      }
      printf("Total: %d points, shift %d\n", totalvert, boxshift);
      fflush(stdout); #--#
    }
  fclose(outfile);


}

#ifndef NOASSERT
  printf("\nDone\n\n");
  fid = fopen(argv[2],"r");
  fgets(line, 132, fid);
  printf("%s\n", line);
  while (fread(&j,sizeof(int), 1, fid) == 1) {
    printf("level has %d\n", j);
    if (j < 1)
      exit(1);
    for (i = 0; i < j; i++) {
      k = fread(&v, sizeof(pointType), 1, fid);
      if (k != 1)
	printf("ERROR on read pt %d, %d\n", k, ferror(fid));
      #-      else
		printf("  %d %x (%d %x %x %x)f <<%d\n", v.boxno, v.boxno,
	       v.index, v.x, v.y, v.z, BBbits(bb)); #--#
    }
  }
  fclose(fid);#--#
#endif
  return 1;
}

const VERT_OUT 1
const EDGE_OUT 2
const TRI_OUT  4
const TET_OUT  8
const INF_OUT 16 # output infinite vertex, edges, ...

void d3output(FILE* ostream, int select,
  ppointType vert, int nvert, pcornerType s, int ncorners) {
int j = 0;
int infin = (0 == select & INF_OUT);

if (0 != select&VERT_OUT) {
if (infin)
     fprintf(ostream, "(%5d:  oo  %5.0f %5.0f %5.0f %10.0f)\n", j,
       vert[j].x, vert[j].y, vert[j].z, vert[j].sq);
  for (i = 1; i < nvert; i++)
     fprintf(ostream, "(%5d: %5d %5.0f %5.0f %5.0f %10.0f)\n", j,
    vert[j].index, vert[j].x, vert[j].y, vert[j].z, vert[j].sq);
     }

if (0 != select&TET_OUT) {
for (j = 0; j < ncorners; j += 4)
     if (infin || (s[j].v > vert && s[j+1].v > vert))
     fprintf(ostream, "%5d %5d %5d %5d\n", s[j].v-vert, s[j+1].v-vert, s[j+2]-vert, s[j+3]-vert);
}
}


end # module Tet3
