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
void doAffineAlignment(AffineAlignObj&,const SimMatrix& s, double go, double ge, bool OverlapAlignment);

void getAffineAlignedIndices(AffineAlignObj &affineAlignObj);

double getOlapAffineAlignStartIndices(double* MatrixM, double* MatrixA, double* MatrixB, int ROW_SIZE, int COL_SIZE, int &OlapStartRow, int &OlapStartCol, tbJump &MatrixName);
} // namespace DIAlign

#endif // AFFINEALIGNMENT_H

