

#include "hcurvedata.h"
#include "stdio.h"

int main(int argc, char* argv[])
{
  /// const int pDim = 2;
  /// const int pDim = 3;
  /// const int pDim = 4;
  const int pDim = 5;
  _InitializeHcurve(pDim);

  printf("_KD_ROTATIONS_SIZE == %d\n", _KD_ROTATIONS_SIZE(pDim));

  int num_orthant = (1u<<pDim);
  for (int ort = 0; ort < num_orthant*2; ort++)
    printf("%c", rotations[ort]);
  printf("\n");

  for (int ort = 0; ort < num_orthant; ort++)
    printf("%d ", HILBERT_TABLE[ort]);
  printf("\n");

  return 0;
}
