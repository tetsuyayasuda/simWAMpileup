//--------------------------------------------
//     readSpecWAM.hh
//--------------------------------------------
// 2013-09-13 yasuda

#ifndef READSPECWAM_HH
#define READSPECWAM_HH

#include <string>
#include <vector>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

class readSpecWAM
{
public:
  readSpecWAM();
  ~readSpecWAM();
  void setSpecFits(const std::string inputName);
  double getMatrixCount(int channelNum);
  int getMatrixColNum();
  int getMatrixRowNum();
  sli::mdarray_double getSpecMatrix();
  double getLivetime();
  double getExposure();

private:
  sli::fitscc fits;
  ssize_t sz;

  sli::mdarray_double matrix;
  int matrixColNum;
  int matrixRowNum;
  double livetime;
  double exposure;
  //std::vector<double> integralSpec;

};

#endif
