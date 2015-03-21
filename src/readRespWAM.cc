//--------------------------------------------
//             baseRespWAM.hh
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

#include "readRespWAM.hh"

#define PI 3.141592653589793

using namespace sli;

///////  Public Functions  ////////
readRespWAM::readRespWAM(){
  matrix.init( false );
}

readRespWAM::~readRespWAM(){
}

void readRespWAM::setRespFits(const std::string inputName){
  sz = fits.read_stream( inputName.c_str() );
  if( sz < 0 ){
    std::cout << "Could not read " << inputName << std::endl;
    exit(1);
  }

  std::cout << "Response successfully loaded: " << inputName << std::endl;  
  matrixColNum = fits.table("MATRIX").col_length();
  eboundsColNum = fits.table("EBOUNDS").col_length();
  matrixRowNum = fits.table("MATRIX").row_length();
  eboundsRowNum = fits.table("EBOUNDS").row_length();

  fits_table_col fitsMatrix = fits.table("MATRIX").col("MATRIX");
  matrix.resize_2d(eboundsRowNum, eboundsRowNum);
  for(int j=0; j<eboundsRowNum; j++){ // X
  
    for(int k=0; k<eboundsRowNum; k++){ // Y
      double inteResp = 0;
      double eboundsLo = getEboundsRangeLow( k );
      double eboundsHi = getEboundsRangeHigh( k );
    
      for (int i=0; i<matrixRowNum; i++) { 	
	double matrixLo = getMatrixRangeLow( i );
	double matrixHi = getMatrixRangeHigh( i );
	
	if( matrixLo <= eboundsLo && eboundsLo < matrixHi){
	  if( eboundsHi <= matrixHi ){
	    inteResp += (eboundsHi - eboundsLo) * fitsMatrix.dvalue(i, j);
	  }
	  else{
	    inteResp += (matrixHi - eboundsLo) * fitsMatrix.dvalue(i, j);
	  }
	}
	else if( eboundsLo <= matrixLo && matrixLo < eboundsHi ){
	  if( matrixHi <= eboundsHi ){
	    inteResp += (matrixHi - matrixLo) * fitsMatrix.dvalue(i, j);
	  }
	  else{
	    inteResp += (eboundsHi - matrixLo) * fitsMatrix.dvalue(i, j);
	  }
	}
	else if( eboundsHi <= matrixLo){
	  break;
	}
      }
      matrix(k, j) = inteResp / (eboundsHi - eboundsLo);
    }
  }
}

int readRespWAM::getMatrixColNum(){return matrixColNum;}
int readRespWAM::getMatrixRowNum(){return matrixRowNum;}
int readRespWAM::getEboundsColNum(){return eboundsColNum;}
int readRespWAM::getEboundsRowNum(){return eboundsRowNum;}

double readRespWAM::getMatrixRangeLow(int channelNum){return fits.table("MATRIX").col("ENERG_LO").dvalue(channelNum);}
double readRespWAM::getMatrixRangeHigh(int channelNum){return fits.table("MATRIX").col("ENERG_HI").dvalue(channelNum);}
double readRespWAM::getEboundsRangeLow(int channelNum){return fits.table("EBOUNDS").col("E_MIN").dvalue(channelNum);}
double readRespWAM::getEboundsRangeHigh(int channelNum){return fits.table("EBOUNDS").col("E_MAX").dvalue(channelNum);}

sli::mdarray_double readRespWAM::getRespMatrix(){return matrix;}

