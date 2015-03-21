//--------------------------------------------
//     readRespWAM.hh
//--------------------------------------------
// 2013-09-13 yasuda

#ifndef READRESPWAM_HH
#define READRESPWAM_HH

#include <string>
#include <vector>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

class readRespWAM
{
public:
  readRespWAM();
  ~readRespWAM();
  void setRespFits(const std::string inputName);
  int getColumnNum();
  double getMatrixRangeLow(int channelNum);
  double getMatrixRangeHigh(int channelNum);
  double getEboundsRangeLow(int channelNum);
  double getEboundsRangeHigh(int channelNum);
  int getMatrixColNum();
  int getMatrixRowNum();
  int getEboundsColNum();
  int getEboundsRowNum();
  sli::mdarray_double getRespMatrix();
  //sli::mdarray_double getRespReverseMatrix();

private:
  sli::fitscc fits;
  ssize_t sz;

  sli::mdarray_double matrix;
  int matrixColNum, eboundsColNum;
  int matrixRowNum, eboundsRowNum;

  
};

#endif
