//--------------------------------------------
//             readSpecWAM.hh
//--------------------------------------------
// 2013-09-13 yasuda

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

#include "readSpecWAM.hh"

#define PI 3.141592653589793

using namespace sli;

///////  Public Functions  ////////
readSpecWAM::readSpecWAM(){
  matrix.init( false );
}

readSpecWAM::~readSpecWAM(){
}

void readSpecWAM::setSpecFits(const std::string inputName){
  sz = fits.read_stream( inputName.c_str() );
  if( sz < 0 ){
    std::cout << "Could not read " << inputName << std::endl;
    exit(1);
  }

  std::cout << "Spectrum successfully loaded: " << inputName << std::endl;  
  matrixColNum = fits.table("SPECTRUM").col_length();
  matrixRowNum = fits.table("SPECTRUM").row_length();
  livetime = fits.hdu("SPECTRUM").headerf("LIVETIME").dvalue();
  exposure = fits.hdu("SPECTRUM").headerf("EXPOSURE").dvalue();
  if( isnan( livetime ) != 0 )livetime = 1.0;
  std::cout << " livetime(sec): " << livetime << std::endl;
  std::cout << " exposure(sec): " << exposure << std::endl;
  
  int rateFlag = 0;
  matrix.resize_2d(matrixRowNum, 1);
  for (int i=0; i<matrixRowNum; i++) {
    matrix(i, 0) = fits.table("SPECTRUM").col("COUNTS").dvalue(i);
    if( !fits.table("SPECTRUM").col("COUNTS").dvalue(i) ){
      matrix(i, 0) = fits.table("SPECTRUM").col("RATE").dvalue(i);
      rateFlag = 1;
    }
  }
  if( rateFlag == 1 ){
    std::cout << " RATE has been recorded. exposure(sec): " << exposure << " -> 1.0" << std::endl;
    exposure = 1.0;
  } 
}

int readSpecWAM::getMatrixColNum(){return matrixColNum;}
int readSpecWAM::getMatrixRowNum(){return matrixRowNum;}

double readSpecWAM::getMatrixCount(int channelNum){return fits.table("SPECTRUM").col("COUNTS").dvalue(channelNum);}

sli::mdarray_double readSpecWAM::getSpecMatrix(){return matrix;}

double readSpecWAM::getLivetime(){return livetime;};
double readSpecWAM::getExposure(){return exposure;};
