/*----------------------------------------------------------------------------
The MIT License (MIT)

  Copyright (c) 2015 Tetsuya Yasuda

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in 
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
---------------------------------------------------------------------------*/

//************************************
//   simWAMpileup
//   T.Yasuda   2013-09-13 ver0.0
//   T.Yasuda   2013-09-25 ver0.1
//************************************

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

#include "readRespWAM.hh"
#include "readSpecWAM.hh"
#include "simPileup.hh"


sli::mdarray_double innerProduct(sli::mdarray_double model, sli::mdarray_double resp, int eBoundsNum);
void makeSpecTemplate(char* outName, int energyBins, double deadTime);


int main(int argc, char *argv[]){

  // rundom seed
  timeval time_val;
  gettimeofday( &time_val, NULL );
  int seed = (unsigned)time(NULL) + time_val.tv_usec;
  gRandom->SetSeed( seed );


  if( argc != 7 ){
    std::cout << "Usage: ./readResponse [resp] [bgd] [src] [incidentCountRate] [lightcurve.qdp] [overFlowRate]" 
	      << std::endl << std::endl;
    return 0;
  }
  
  // initialize  ------------------------------------
  char* respName = argv[1];
  char* bgdName = argv[2];
  char* srcName = argv[3];
  int incidentCount = atoi( argv[4] );
  char* lcName = NULL;
  if( argv[5] != NULL ) lcName = argv[5];
  int overflowRate = 0;
  if( argv[6] != NULL ) overflowRate = atoi( argv[6] );
  clock_t startTime, endTime;
  startTime = clock();
  

  // read response matrix ---------------------------
  readRespWAM* readResp = new readRespWAM();
  readResp->setRespFits( respName );
  int energyBins = readResp->getEboundsRowNum();
  std::cout << "readResp ........... done" << std::endl;


  // read bgd fits ----------------------------------
  readSpecWAM* readBgd = new readSpecWAM();
  readBgd->setSpecFits( bgdName );
  double exposureBgd = readBgd->getExposure();
  sli::mdarray_double bgdMatrix = readBgd->getSpecMatrix();
  bgdMatrix /= exposureBgd;
  std::cout << "readBgd ........... done" << std::endl;


  // read src fits (fakeit pha) --------------
  readSpecWAM* readSrc = new readSpecWAM();
  readSrc->setSpecFits( srcName );
  double exposureSrc = readSrc->getExposure();
  sli::mdarray_double srcMatrix = readSrc->getSpecMatrix();
  srcMatrix /= exposureSrc;
  std::cout << "readSrc ........... done" << std::endl;


  // simPileup --------------------------------------
  simPileup* runPileup = new simPileup();
  runPileup->setCountRate( incidentCount );
  for(int i=0; i<energyBins; i++){
    runPileup->setEnergyRange( readResp->getEboundsRangeLow( i ) );
  }
  runPileup->setEnergyRange( readResp->getEboundsRangeHigh( energyBins-1 ));
  runPileup->preparation();
  if(lcName != NULL){
    runPileup->inteLightcurve( lcName );
    std::cout << "inteLightcurve ........... done (" << lcName << ")" <<std::endl;}
  else{
    runPileup->inteLightcurve( );
    std::cout << "inteLightcurve ........... done" <<std::endl;}
  runPileup->inteBgd( bgdMatrix );
  runPileup->inteSrc( srcMatrix );
  runPileup->setOverflowRate( overflowRate, 1.0 );
  runPileup->runPileup();
  std::cout << "runPileup ......... done" << std::endl;
  double deadTime = runPileup->getDeadtime();
  double liveRate = runPileup->getLiveCount();


  // get matricies  ---------------------------------- 
  sli::mdarray_double randOriginMatrix = runPileup->getOriginMatrix();
  sli::mdarray_double randSrcMatrix = runPileup->getSrcMatrix();
  sli::mdarray_double randBgdMatrix = runPileup->getBgdMatrix();
  sli::mdarray_double randDeadMatrix = runPileup->getDeadMatrix();


  // make FITS header  ---------------------------------- 
  std::string templateOrigin = "Origin";
  std::string templateSrc = "Src";
  std::string templateBgd = "Bgd";
  std::string templateDead = "Dead";
  makeSpecTemplate( (char*)("template"+templateOrigin+".txt").c_str(), energyBins, 0.0);
  makeSpecTemplate( (char*)("template"+templateSrc+".txt").c_str(), energyBins, 0.0);
  makeSpecTemplate( (char*)("template"+templateBgd+".txt").c_str(), energyBins, 0.0);
  makeSpecTemplate( (char*)("template"+templateDead+".txt").c_str(), energyBins, deadTime);
  

  // define fits ----------------------------------------
  sli::fitscc output0, output1, output2, output3;
  ssize_t sz0, sz1, sz2, sz3;
  output0.read_template( (char*)("template"+templateOrigin+".txt").c_str() );
  output1.read_template( (char*)("template"+templateSrc+".txt").c_str() );
  output2.read_template( (char*)("template"+templateBgd+".txt").c_str() );
  output3.read_template( (char*)("template"+templateDead+".txt").c_str() );


  // write spectrum fits ------------------------------
  for(int i=0; i<energyBins; i++){
    output0.table("SPECTRUM").col("CHANNEL").assign( i+1, i );
    output0.table("SPECTRUM").col("COUNTS").assign( randOriginMatrix.dvalue(i, 0), i );
    output1.table("SPECTRUM").col("CHANNEL").assign( i+1, i );
    output1.table("SPECTRUM").col("COUNTS").assign( randSrcMatrix.dvalue(i, 0), i );
    output2.table("SPECTRUM").col("CHANNEL").assign( i+1, i );
    output2.table("SPECTRUM").col("COUNTS").assign( randBgdMatrix.dvalue(i, 0), i );
    output3.table("SPECTRUM").col("CHANNEL").assign( i+1, i );
    output3.table("SPECTRUM").col("COUNTS").assign( randDeadMatrix.dvalue(i, 0), i );
  }
  sz0 = output0.write_stream( (char*)(templateOrigin+".fits").c_str() );  
  sz1 = output1.write_stream( (char*)(templateSrc+".fits").c_str() );  
  sz2 = output2.write_stream( (char*)(templateBgd+".fits").c_str() );  
  sz3 = output3.write_stream( (char*)(templateDead+".fits").c_str() );  
  std::cout << "FITS files ....... created" << std::endl;


  // output simulation result --------------------------
  endTime = clock();
  time_t time_now;
  time( &time_now );
  struct tm* tm_now =localtime( &time_now );
  std::ofstream output_result;
  output_result.open("output_result.dat");
  output_result << "#########  Result  #########" << std::endl;
  output_result << "Date:              " << tm_now->tm_year+1900 <<"-"<< tm_now->tm_mon+1 <<"-"<< tm_now->tm_mday <<"_"
		<< tm_now->tm_hour<<":"<< tm_now->tm_min <<":"<< tm_now->tm_sec << std::endl;
  output_result << "ElapsedTime[sec]:  "<<  (double)(endTime - startTime)/CLOCKS_PER_SEC << std::endl;
  output_result << "Random seed:       "<< seed << std::endl;
  output_result << "---------  Input  ----------" << std::endl;
  output_result << "Src.Fits(fakeit):  "<< srcName << std::endl;
  output_result << "Bgd.Fits:          "<< bgdName << std::endl;
  output_result << "Resp.Fits:         "<< respName << std::endl;
  if( lcName != NULL ) output_result << "LightCurve: "<< lcName << std::endl;
  output_result << "IncidentRate[c/s]: "<< incidentCount << std::endl;
  output_result << "Exposure[sec]:     "<< 1.0 << std::endl;
  output_result << "---------  Output  ---------" << std::endl;
  output_result << "Src.Pha:           " << templateSrc+".fits" << std::endl;
  output_result << "Bgd.Pha:           " << templateBgd+".fits" << std::endl;
  output_result << "Src+Bgd.Pha:       " << templateOrigin+".fits" << std::endl;
  output_result << "Pileuped.Pha:      " << templateDead+".fits" << std::endl;
  output_result << "DeadTime[sec]:     " << deadTime << std::endl;
  output_result << "LiveTime[sec]:     " << 1.0-deadTime << std::endl;
  output_result << "LiveRate[c/s]:     " << liveRate << std::endl;
  output_result << "NormDeadCorreFact: " << 1.0/(1.0-deadTime) << std::endl;
  output_result << "PileupCorrectFact: " << (double)incidentCount/(double)liveRate << std::endl;
  output_result.close();



  delete readResp;
  delete readBgd;
  delete readSrc;
  delete runPileup;

  return 0;
}

sli::mdarray_double innerProduct(sli::mdarray_double model, sli::mdarray_double resp, int eBoundsNum){
  sli::mdarray_double result;
  result.init(false);
  result.resize_2d(eBoundsNum, 1);
  for(int j=0; j<eBoundsNum; j++){
    result(j, 0) = 0;
    for(int i=0; i<eBoundsNum; i++){
      result(j, 0) += model.dvalue(j, 0) * resp.dvalue(j, i);
      std::cout << resp.dvalue(j, i) <<" ";
    }
    std::cout << std::endl;
  }
  return result;
}


void makeSpecTemplate(char* outName, int energybins, double deadTime){
  int binLo, binHi;
  if(energybins == 54){
    binLo = 1;
    binHi = 54;
  }
  else{
    binLo = 0;
    binHi = 54;
  }
  
  std::ofstream output;
  output.open( outName );
  output << "SIMPLE  =                    T / file does conform to FITS standard " << std::endl;
  output << "BITPIX  =                  -32 / number of bits per data pixel " << std::endl;
  output << "NAXIS   =                    0 / number of data axes " << std::endl;
  output << "EXTEND  =                    T / FITS dataset may contain extensions " << std::endl;
  output << "COMMENT   FITS (Flexible Image Transport System) format is defined in \'Astronomy " << std::endl;
  output << "COMMENT   and Astrophysics\', volume 376, page 359; bibcode: 2001A&A...376..359H " << std::endl;
  output << "END " << std::endl;
  output << std::endl;
  output << "XTENSION= \'BINTABLE\'           / binary table extension " << std::endl;
  output << "BITPIX  =                    8 / 8-bit bytes " << std::endl;
  output << "NAXIS   =                    2 / 2-dimensional binary table " << std::endl;
  output << "NAXIS1  =                   10 / width of table in bytes " << std::endl;
  output << "NAXIS2  =                   "<< binHi <<" / number of rows in table " << std::endl;
  output << "PCOUNT  =                    0 / size of special data area " << std::endl;
  output << "GCOUNT  =                    1 / one data group (required keyword) " << std::endl;
  output << "TFIELDS =                    3 / number of fields in each row " << std::endl;
  output << "TTYPE1  = \'CHANNEL \'           / Detector channel (type unknown) " << std::endl;
  output << "TFORM1  = \'J       \'           / data format of field: 4-byte INTEGER " << std::endl;
  output << "TTYPE2  = \'COUNTS  \'           / Counts per channel " << std::endl;
  output << "TFORM2  = \'E       \'           / data format of field: 4-byte REAL " << std::endl;
  output << "TUNIT2  = \'count   \'           / physical unit of field " << std::endl;
  output << "TTYPE3  = \'QUALITY \'           / Quality flag of this channel (0=good) " << std::endl;
  output << "TFORM3  = \'I       \'           / data format of field: 2-byte INTEGER " << std::endl;
  output << "EXTNAME = \'SPECTRUM\'           / name of this binary table extension " << std::endl;
  output << "HDUCLASS= \'OGIP    \'           / format conforms to OGIP standard " << std::endl;
  output << "HDUCLAS1= \'SPECTRUM\'           / PHA dataset (OGIP memo OGIP-92-007) " << std::endl;
  output << "HDUVERS1= \'1.2.0   \'           / Obsolete - included for backwards compatibility" << std::endl;
  output << "HDUVERS = \'1.2.0   \'           / Version of format (OGIP memo OGIP-92-007)" << std::endl;
  output << "HDUCLAS2= \'DERIVED \'           / WARNING This is NOT an OGIP-approved value" << std::endl;
  output << "HDUCLAS3= \'COUNT   \'           / PHA data stored as Counts (not count/s)" << std::endl;
  output << "TLMIN1  =                    "<< binLo <<" / Lowest legal channel number" << std::endl;
  output << "TLMAX1  =                   "<< binHi <<" / Highest legal channel number" << std::endl;
  output << "TELESCOP= \'SUZAKU  \'           / mission/satellite name" << std::endl;
  output << "INSTRUME= \'HXD     \'           / instrument/detector name" << std::endl;
  output << "DETNAM  = \'WAM_ANTI\'           / specific detector name in use" << std::endl;
  output << "FILTER  = \'NONE    \'           / filter in use" << std::endl;
  output << "EXPOSURE=         "<< 1.0-deadTime <<" / exposure (in seconds)" << std::endl;
  output << "AREASCAL=         1.000000E+00 / area scaling factor" << std::endl;
  output << "BACKFILE= \'NONE    \'           / associated background filename" << std::endl;
  output << "BACKSCAL=         1.000000E+00 / background file scaling factor" << std::endl;
  output << "CORRFILE= \'NONE    \'           / associated correction filename" << std::endl;
  output << "CORRSCAL=         1.000000E+00 / correction file scaling factor" << std::endl;
  output << "RESPFILE= \'NONE    \'           / associated redistrib matrix filename" << std::endl;
  output << "ANCRFILE= \'NONE    \'           / associated ancillary response filename" << std::endl;
  output << "PHAVERSN= \'1992a   \'           / obsolete" << std::endl;
  output << "DETCHANS=                   "<< binHi <<" / total number possible channels" << std::endl;
  output << "CHANTYPE= \'UNKNOWN \'           / channel type (PHA, PI etc)" << std::endl;
  output << "POISSERR=                    T / Poissonian errors to be assumed" << std::endl;
  output << "STAT_ERR=                    0 / no statistical error specified" << std::endl;
  output << "SYS_ERR =                    0 / no systematic error specified" << std::endl;
  output << "GROUPING=                    0 / no grouping of the data has been defined" << std::endl;
  output << "END" << std::endl;
  output.close();
}
