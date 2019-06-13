#ifndef AFFINEALIGNMENT_H
#define AFFINEALIGNMENT_H

#include "affinealignobj.h"
#include "alignment.h"
#include "utils.h"
#include "similarityMatrix.h"
#include <limits>

namespace DIAlign 
{
// It performs affine alignment on similarity matrix and fills three matrices M, A and B, and corresponding traceback matrices.
AffineAlignObj doAffineAlignment(SimMatrix s, double go, double ge, bool OverlapAlignment);

void getAffineAlignedIndices(AffineAlignObj &affineAlignObj);

template<class T>
double getOlapAffineAlignStartIndices(T MatrixM, T MatrixA, T MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);
} // namespace DIAlign

#endif // AFFINEALIGNMENT_H

