//-----------------------
//   simPileup
//
//   Tetsuya Yasuda
//-----------------------

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TObjString.h>
#include <TH1D.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TLinearMinimizer.h>
#include <Math/WrappedTF1.h>

#include "readText.hh"
#include "powCalc.hh"
#include "simPileup.hh"

//using namespace std;

simPileup::simPileup(){
  // parameters ------------------
  nearEventNum = 1000;
  numBounds = 0;
  counts = 0;
  lcBins = 0;
  specBins = 0;
  exposure = 0.0;
  countsBgd = 0;
  countDeadTime = 0;
  countLiveNumber = 0;
  separateStepLc = 30;
  separateStepSpec = 20;
  overflowCounts = 0;
}


simPileup::~simPileup(){
}


void simPileup::setEnergyRange(double eBounds){
  v_energyBounds.push_back( eBounds );
  if( numBounds == 0 ){
    v_energyBoundsQdp.push_back( eBounds );
  }
  else{
    for(int j=0; j<separateStepSpec; j++){
      double prePoint = v_energyBounds[ numBounds-1 ];
      double postPoint = (eBounds - prePoint)/((double)separateStepSpec)*(j+1) +prePoint;
      v_energyBoundsQdp.push_back( postPoint );
    }
  }
  numBounds++;
  threshEnergy = v_energyBounds[0];
  threshEneUpp = v_energyBounds[ v_energyBounds.size()-1 ];
}


void simPileup::setCountRate(int countRate){
  counts = countRate;
}


void simPileup::preparation(){
  if( counts == 0 && numBounds < 2 ){
    std::cout << "ERROR: COUNTS or NUMBOUNDS have not been set." << std::endl;
    exit(1);
  }
  hCount = new TH1D("Count", "Spectrum", v_energyBounds.size()-1, (double*)&(v_energyBounds.at(0)));
  hCountBgd = new TH1D("CountBgd", "Spectrum", v_energyBounds.size()-1, (double*)&(v_energyBounds.at(0)));
  hCountSrc = new TH1D("CountSrc", "Spectrum", v_energyBounds.size()-1, (double*)&(v_energyBounds.at(0)));
  hCountDead = new TH1D("CountDead", "Spectrum", v_energyBounds.size()-1, (double*)&(v_energyBounds.at(0)));
  
  hCountQdp = new TH1D("CountQdp", "Spectrum", v_energyBoundsQdp.size()-1, (double*)&(v_energyBoundsQdp.at(0)));
  hCountBgdQdp = new TH1D("CountBgdQdp", "Spectrum", v_energyBoundsQdp.size()-1, (double*)&(v_energyBoundsQdp.at(0)));
  hCountSrcQdp = new TH1D("CountSrcQdp", "Spectrum", v_energyBoundsQdp.size()-1, (double*)&(v_energyBoundsQdp.at(0)));
  hCountDeadQdp = new TH1D("CountDeadQdp", "Spectrum", v_energyBoundsQdp.size()-1, (double*)&(v_energyBoundsQdp.at(0)));
}


void simPileup::inteLightcurve(){
  v_inteLc.resize( 2 );
  v_inteStepLc.resize( 2 );
  v_inteLc[0] = 0.0;  
  v_inteLc[1] = 1.0;  
  v_inteStepLc[0] = 0.0;  
  v_inteStepLc[1] = 1.0;  
  lcBins = v_inteLc.size();
  exposure = 1.0;
  hLightCurve = new TH1D("LightCurve", "Lightcurve", 200, 0, exposure);
}


void simPileup::inteLightcurve(char* lcName){
  readText* lc = new readText();
  lc->setFile( lcName );
  lc->setColumn2();
  int lcEventNum = lc->getEventNum();

  // integral (light curve) -----------------------------------
  v_inteLc.push_back( 0.0 );
  v_inteStepLc.push_back( lc->getColumn1(0) );
  double inteSumLc = 0.0;
  for(int i=0; i<lcEventNum-1; i=i+1){
    double X1 = lc->getColumn1(i) - lc->getColumn1(0);
    double Y1 = lc->getColumn2(i);
    double X2 = lc->getColumn1(i+1) - lc->getColumn1(0);
    double Y2 = lc->getColumn2(i+1);
    
    if(Y1 < 10e-256 || Y2 < 10e-256){
      std::cout << X1 <<" "<< X2 <<" "<< lcEventNum << std::endl;
      std::cout << "ERROR: "<< Y1 <<" or "<< Y2 << " is too small to calculate." << std::endl;
      exit(1);
    }
    double A = linearConst(X1, Y1, X2, Y2);
    double B = linearIndex(X1, Y1, X2, Y2);
    for(int j=0; j<separateStepLc; j++){
      double prePoint = (X2-X1)/((double)separateStepLc)*(j) + X1;
      double postPoint = (X2-X1)/((double)separateStepLc)*(j+1) + X1;
      double integral = linearInte(prePoint, postPoint, A, B);
      inteSumLc = inteSumLc + integral;
      v_inteLc.push_back( inteSumLc );
      v_inteStepLc.push_back( postPoint );
    }
  }
  lcBins = v_inteLc.size();
  exposure = v_inteStepLc[ lcBins-1 ] - v_inteStepLc[ 0 ];
  /*
  ofstream output_inteSrcQdp;
  output_inteSrcQdp.open("output_inteLcSrc.qdp");
  output_inteSrcQdp << "line on" << std::endl;
  output_inteSrcQdp << "lw 5" << std::endl;
  output_inteSrcQdp << "la x Time (sec)" << std::endl;
  for(int i=0; i<lcBins; i++){
    output_inteSrcQdp << v_inteStepLc[i] <<" "<< v_inteLc[i] << std::endl;
  }
  output_inteSrcQdp.close(); 
  std::cout << "output_inteSpecSrc.qdp has been created." <<std::endl;
  */
  hLightCurve = new TH1D("LightCurve", "Lightcurve", 200, 0, exposure);
}

 
void simPileup::inteSrc(sli::mdarray_double matrix){
  matrixSrc.init( false );
  matrixSrc = matrix;
  int specEventNum = matrixSrc.length();
  
  // integral (spectrum) -----------------------------------
  double inteSumSrc = 0.0;
  int skipFlag = 0;
  v_inteSrc.push_back( inteSumSrc );
  v_inteStepSrc.push_back( v_energyBounds[0] );
  //inteSumSrc += (0.5*(v_energyBounds[1]-v_energyBounds[0])) * matrixSrc.dvalue(0);
  //v_inteSrc.push_back( inteSumSrc );
  //v_inteStepSrc.push_back( (0.5*(v_energyBounds[1]-v_energyBounds[0]))+v_energyBounds[0] );

  for(int i=0; i<specEventNum-1; i++){
    double X1 = ((v_energyBounds[ i+1 ]-v_energyBounds[ i ]) *0.5) + v_energyBounds[ i ];
    double X2 = ((v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]) *0.5) + v_energyBounds[ i+1 ];
    double Y1 = matrixSrc.dvalue( i )/(v_energyBounds[ i+1 ]-v_energyBounds[ i ]);
    double Y2 = matrixSrc.dvalue( i+1 )/(v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]);
    
    if(Y1 < 10e-256 || Y2 < 10e-256){
      cout <<"ERROR: "<< Y1 <<" or "<< Y2 << " is too small to calculate." << endl;
      skipFlag = 1;
    }

    double A = powConst(X1, Y1, X2, Y2);
    double B = powIndex(X1, Y1, X2, Y2);

    if(i == 0){
      double X00 = v_energyBounds[ 0 ];
      double X01 = X1;
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = powInte(prePoint, postPoint, A, B);
	if( isnan( integral ) || isinf( integral ) ) skipFlag = 1; 
	if( skipFlag == 1 ) integral = 0.0;
	inteSumSrc += integral;
	v_inteSrc.push_back( inteSumSrc );
	v_inteStepSrc.push_back( postPoint );
      }
    }
    
    for(int j=0; j<separateStepSpec; j++){
      double prePoint = (X2-X1)/((double)separateStepSpec)*(j) + X1;
      double postPoint = (X2-X1)/((double)separateStepSpec)*(j+1) + X1;
      double integral = powInte(prePoint, postPoint, A, B);
      if( isnan( integral ) || isinf( integral ) ) skipFlag = 1; 
      if( skipFlag == 1 ) integral = 0.0;
      inteSumSrc += integral;
      v_inteSrc.push_back( inteSumSrc );
      v_inteStepSrc.push_back( postPoint );
    }

    if(i == specEventNum-2){
      double X00 = X2;
      double X01 = v_energyBounds[ specEventNum ];
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = powInte(prePoint, postPoint, A, B);
	if( isnan( integral ) || isinf( integral ) ) skipFlag = 1; 
	if( skipFlag == 1 ) integral = 0.0;
	inteSumSrc += integral;
	v_inteSrc.push_back( inteSumSrc );
	v_inteStepSrc.push_back( postPoint );
      }
    }
    skipFlag = 0;
  }

  //inteSumSrc += (0.5*(v_energyBounds[specEventNum]-v_energyBounds[specEventNum-1])) * matrixSrc.dvalue(specEventNum-1);
  //v_inteSrc.push_back( inteSumSrc );
  //v_inteStepSrc.push_back( v_energyBounds[specEventNum] );
  
  specBins = v_inteSrc.size();
  /*
  ofstream output_inteSrcQdp;
  output_inteSrcQdp.open("output_inteSpecSrc.qdp");
  output_inteSrcQdp << "log x" << std::endl;
  output_inteSrcQdp << "line on" << std::endl;
  output_inteSrcQdp << "lw 5" << std::endl;
  output_inteSrcQdp << "la x Energy (keV)" << std::endl;
  for(int i=0; i<specBins; i++){
    output_inteSrcQdp << v_inteStepSrc[i] <<" "<< v_inteSrc[i] << std::endl;
  }
  output_inteSrcQdp.close(); 
  std::cout << "output_inteSpecSrc.qdp has been created." << std::endl;
  */
}


void simPileup::inteBgd(sli::mdarray_double matrix){
  matrixBgd.init( false );
  matrixBgd = matrix;
  int specEventNum = matrixBgd.length();
  
  // integral (spectrum) -----------------------------------
  v_inteBgd.push_back( 0.0 );
  v_inteStepBgd.push_back( v_energyBounds[0] );
  double inteSumBgd = 0.0;
  double preCountsBgd = 0.0;
  int skipFlag = 0;
  //inteSumBgd += (0.5*(v_energyBounds[1]-v_energyBounds[0])) * matrixBgd.dvalue(0);
  //v_inteBgd[1] = inteSumBgd;
  //v_inteStepBgd[1] = (0.5*(v_energyBounds[1]-v_energyBounds[0]))+v_energyBounds[0];

  for(int i=0; i<specEventNum-1; i=i+1){
    double X1 = ((v_energyBounds[ i+1 ]-v_energyBounds[ i ]) *0.5) + v_energyBounds[ i ];
    double X2 = ((v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]) *0.5) + v_energyBounds[ i+1 ];
    double Y1 = matrixBgd.dvalue( i )/(v_energyBounds[ i+1 ]-v_energyBounds[ i ]);
    double Y2 = matrixBgd.dvalue( i+1 )/(v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]);
    
    if(Y1 < 10e-256 || Y2 < 10e-256){
      //cout <<"ERROR: "<< Y1 <<"or"<< Y2 << " is too small to calculate." << endl;
      skipFlag = 1;
    }
    
    double A = linearConst(X1, Y1, X2, Y2);
    double B = linearIndex(X1, Y1, X2, Y2);

    if(i == 0){
      double X00 = v_energyBounds[ 0 ];
      double X01 = X1;
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = linearInte(prePoint, postPoint, A, B);
	if( skipFlag == 1 ) integral = 0.0;
	inteSumBgd += integral;
	v_inteBgd.push_back( inteSumBgd );
	v_inteStepBgd.push_back( postPoint );
      }
    }

    for(int j=0; j<separateStepSpec; j++){
      double prePoint = (X2-X1)/((double)separateStepSpec)*(j) + X1;
      double postPoint = (X2-X1)/((double)separateStepSpec)*(j+1) + X1;
      double integral = linearInte(prePoint, postPoint, A, B);
      if( skipFlag == 1 ) integral = 0.0;
      inteSumBgd += integral;
      v_inteBgd.push_back( inteSumBgd );
      v_inteStepBgd.push_back( postPoint );
    }

    if(i == specEventNum-2){
      double X00 = X2;
      double X01 = v_energyBounds[ specEventNum ];
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = linearInte(prePoint, postPoint, A, B);
	if( skipFlag == 1 ) integral = 0.0;
	inteSumBgd += integral;
	v_inteBgd.push_back( inteSumBgd );
	v_inteStepBgd.push_back( postPoint );
      }
    }
    skipFlag = 0;
    preCountsBgd += matrixBgd.dvalue( i );
  }

  /*
  for(int i=0; i<specEventNum-1; i=i+1){
    double X1 = ((v_energyBounds[ i+1 ]-v_energyBounds[ i ]) *0.5) + v_energyBounds[ i ];
    double X2 = ((v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]) *0.5) + v_energyBounds[ i+1 ];
    double Y1 = matrixBgd.dvalue( i )/(v_energyBounds[ i+1 ]-v_energyBounds[ i ]);
    double Y2 = matrixBgd.dvalue( i+1 )/(v_energyBounds[ i+2 ]-v_energyBounds[ i+1 ]);
    
    if(Y1 < 10e-256 || Y2 < 10e-256){
      //cout <<"ERROR: "<< Y1 <<"or"<< Y2 << " is too small to calculate." << endl;
      skipFlag = 1;
    }
    
    double A = powConst(X1, Y1, X2, Y2);
    double B = powIndex(X1, Y1, X2, Y2);

    if(i == 0){
      double X00 = v_energyBounds[ 0 ];
      double X01 = X1;
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = powInte(prePoint, postPoint, A, B);
	if( skipFlag == 1 ) integral = 0.0;
	inteSumBgd += integral;
	v_inteBgd.push_back( inteSumBgd );
	v_inteStepBgd.push_back( postPoint );
      }
    }

    for(int j=0; j<separateStepSpec; j++){
      double prePoint = (X2-X1)/((double)separateStepSpec)*(j) + X1;
      double postPoint = (X2-X1)/((double)separateStepSpec)*(j+1) + X1;
      double integral = powInte(prePoint, postPoint, A, B);
      if( skipFlag == 1 ) integral = 0.0;
      inteSumBgd += integral;
      v_inteBgd.push_back( inteSumBgd );
      v_inteStepBgd.push_back( postPoint );
    }

    if(i == specEventNum-2){
      double X00 = X2;
      double X01 = v_energyBounds[ specEventNum ];
      for(int j=0; j<separateStepSpec; j++){
	double prePoint = (X01-X00)/((double)separateStepSpec)*(j) + X00;
	double postPoint = (X01-X00)/((double)separateStepSpec)*(j+1) + X00;
	double integral = powInte(prePoint, postPoint, A, B);
	if( skipFlag == 1 ) integral = 0.0;
	inteSumBgd += integral;
	v_inteBgd.push_back( inteSumBgd );
	v_inteStepBgd.push_back( postPoint );
      }
    }
    skipFlag = 0;
    preCountsBgd += matrixBgd.dvalue( i );
  }
  */

  //inteSumBgd += (0.5*(v_energyBounds[specEventNum]-v_energyBounds[specEventNum-1])) * matrixBgd.dvalue(specEventNum-1);
  //v_inteBgd.push_back( inteSumBgd );
  //v_inteStepBgd.push_back( v_energyBounds[specEventNum] );
  preCountsBgd += matrixBgd.dvalue(specEventNum-1);  

  specBins = v_inteBgd.size();
  countsBgd = (int)preCountsBgd;
  //countsBgd = (int)v_inteBgd[ v_inteBgd.size()-1 ];
  /*
  ofstream output_inteBgdQdp;
  output_inteBgdQdp.open("output_inteSpecBgd.qdp");
  output_inteBgdQdp << "log x" << std::endl;
  output_inteBgdQdp << "line on" << std::endl;
  output_inteBgdQdp << "lw 5" << std::endl;
  output_inteBgdQdp << "la x Energy (keV)" << std::endl;
  for(int i=0; i<specBins; i++){
    output_inteBgdQdp << v_inteStepBgd[i] <<" "<< v_inteBgd[i] << std::endl;
  }
  output_inteBgdQdp.close(); 
  */
  std::cout << "output_inteSpecBgd.qdp has been created." << std::endl;
  std::cout << "BGD countrate(c/s): " << countsBgd << std::endl;
}


void simPileup::setOverflowRate(int overflowRate, double liveTimeBgd){
  //if( overflowRate == 0 ) overflowRate = -29.26 + (0.03132 * countsBgd * liveTimeBgd); // countBgd <-- corrected by dead time
  if( overflowRate == 0 ) overflowRate = 350;
  if( overflowRate < 0 ){
    std::cout << "ERROR: OverflowRate is minus. [" << overflowRate <<"]" << std::endl;
    exit(1);
  }
  overflowCounts = (int)overflowRate;
  double overflowStep = 1.0/((double)overflowRate+1.0);
  // exposure = 1.0 sec
  for(int i=1; i<overflowRate+1; i++){
    double holdTime = overflowStep * i;
    v_TandE.push_back(std::pair<double, double>(holdTime, threshEneUpp*3.0));
  }
}


void simPileup::runPileup(){
  if( lcBins == 0 || specBins == 0 ){
    std::cout << "ERROR: LC and/or SPEC have not been set." << std::endl;
    exit(1);
  }
  // count rate ------------------------------
  countsBgd = (int)(exposure * countsBgd);
  countsSrc = counts - countsBgd;
  if( counts <= countsBgd ){
    std::cout << "ERROR: Incident counts is too less than BGD counts." << std::endl;
    exit(1);
  }
  
  // randamize BGD (time, energy) -----------------------------------
  for(int i=0; i<countsBgd; i++){
    double bgdTime = gRandom->Rndm() * exposure;
    double bgdRandInte = gRandom->Rndm() * v_inteBgd[ specBins-1 ];
    
    for(int j=0; j<specBins; j++){
      double Y1 = v_inteBgd[ j ];
      double Y2 = v_inteBgd[ j+1 ];
      if(Y1 < bgdRandInte && bgdRandInte <= Y2){
	double X1 = v_inteStepBgd[ j ];
	double X2 = v_inteStepBgd[ j+1 ];
	double A = linearConst(X1, Y1, X2, Y2);
	double B = linearIndex(X1, Y1, X2, Y2);
	double v = linearInvValue(bgdRandInte, A, B);
	//double A = powConst(X1, Y1, X2, Y2);
	//double B = powIndex(X1, Y1, X2, Y2);
	//double v = powInvValue(bgdRandInte, A, B);
	if( !isnan(v) || !isinf(v) || v<0.1 ){
	  A = linearConst(X1, Y1, X2, Y2);
	  B = linearIndex(X1, Y1, X2, Y2);
	  v = linearInvValue(bgdRandInte, A, B);// [keV]
	}
	v_TandE.push_back(std::pair<double, double>(bgdTime, v));
	hCountBgd->Fill( v, 1.0/exposure );
	hCount->Fill( v, 1.0/exposure );
	hCountBgdQdp->Fill( v, 1.0/exposure );
	hCountQdp->Fill( v, 1.0/exposure );
	hLightCurve->Fill( bgdTime );
	break;
      }
    }
  }
  std::cout << "Randomize BGD (time, energy) ---- done" << std::endl;


  // randamize SRC (time, energy) -----------------------------------
  double srcRandTimeBound = v_inteLc[ lcBins-1 ];
  double srcRandInteBound = v_inteSrc[ specBins-1 ];
  std::cout << "Bin Number Of LightCurve: " << lcBins << std::endl;
  std::cout << "Bin Number Of Spectrum: " << specBins << std::endl;
  int searchBinsWidth = 10;
  int lcBinsStep = lcBins / searchBinsWidth;
  int specBinsStep = specBins / searchBinsWidth;

  for(int i=0; i<countsSrc; i++){
    if( (i+1)%100000 == 0 ){ std::cout << " ( " << i << "/" << countsSrc << " )" << std::endl; }
    double srcRandTime = gRandom->Rndm() * srcRandTimeBound;
    double srcRandInte = gRandom->Rndm() * srcRandInteBound;
    double holdTime, holdEnergy;
    
    int lcBinsFlag = 0;
    for( int j=lcBinsStep; j<lcBins; j+=lcBinsStep ){
      if( srcRandTime <= v_inteLc[ j ] ){
	for( int k=j-lcBinsStep; k<lcBins; k++ ){
	  double Y1 = v_inteLc[ k ];
	  double Y2 = v_inteLc[ k+1 ];
	  if( Y1 < srcRandTime && srcRandTime <= Y2 ){
	    double X1 = v_inteStepLc[ k ];
	    double X2 = v_inteStepLc[ k+1 ];
	    double A = linearConst(X1, Y1, X2, Y2);
	    double B = linearIndex(X1, Y1, X2, Y2);
	    double v = linearInvValue(srcRandTime, A, B);
	    holdTime = v;
	    lcBinsFlag = 1;
	    break;
	  }
	}
      }
      if( lcBinsFlag != 0 ){ break; }
    }
    if( lcBinsFlag == 0 ){
      for(int j=lcBinsStep*searchBinsWidth; j<lcBins; j++){
	double Y1 = v_inteLc[ j ];
	double Y2 = v_inteLc[ j+1 ];
	if( Y1 < srcRandTime && srcRandTime <= Y2 ){
	  double X1 = v_inteStepLc[ j ];
	  double X2 = v_inteStepLc[ j+1 ];
	  double A = linearConst(X1, Y1, X2, Y2);
	  double B = linearIndex(X1, Y1, X2, Y2);
	  double v = linearInvValue(srcRandTime, A, B);
	  holdTime = v;
	  break;
	}
      }
    }
        
    int specBinsFlag = 0;
    for( int j=specBinsStep; j<specBins; j+=specBinsStep ){
      if( srcRandInte <= v_inteSrc[ j ] ){
	for( int k=j-specBinsStep; k<specBins; k++ ){
	  double Y1 = v_inteSrc[ k ];
	  double Y2 = v_inteSrc[ k+1 ];
	  if( Y1 < srcRandInte && srcRandInte <= Y2 ){
	    double X1 = v_inteStepSrc[ k ];
	    double X2 = v_inteStepSrc[ k+1 ];
	    double A = powConst(X1, Y1, X2, Y2);
	    double B = powIndex(X1, Y1, X2, Y2);
	    double v = powInvValue(srcRandInte, A, B);
	    if( !isnan(v) || !isinf(v) || v<0.1 ){
	      A = linearConst(X1, Y1, X2, Y2);
	      B = linearIndex(X1, Y1, X2, Y2);
	      v = linearInvValue(srcRandInte, A, B);
	    }
	    holdEnergy = v;
	    specBinsFlag = 1;
	    break;
	  }
	}
      }
      if( specBinsFlag != 0 ){ break; }
    }
    if( specBinsFlag == 0 ){
      for(int j=specBinsStep*searchBinsWidth; j<specBins; j++){
	double Y1 = v_inteSrc[ j ];
	double Y2 = v_inteSrc[ j+1 ];
	if( Y1 < srcRandInte && srcRandInte <= Y2 ){
	  double X1 = v_inteStepSrc[ j ];
	  double X2 = v_inteStepSrc[ j+1 ];
	  double A = linearConst(X1, Y1, X2, Y2);
	  double B = linearIndex(X1, Y1, X2, Y2);
	  double v = linearInvValue(srcRandInte, A, B);
	  if( !isnan(v) || !isinf(v) || v<0.1 ){
	    A = linearConst(X1, Y1, X2, Y2);
	    B = linearIndex(X1, Y1, X2, Y2);
	    v = linearInvValue(srcRandInte, A, B);
	  }
	  holdEnergy = v;
	  break;
	}
      }
    }
    
    v_TandE.push_back(std::pair<double, double>(holdTime, holdEnergy));
    hCountSrc->Fill( holdEnergy, 1.0/exposure );
    hCount->Fill( holdEnergy, 1.0/exposure );
    hCountSrcQdp->Fill( holdEnergy, 1.0/exposure );
    hCountQdp->Fill( holdEnergy, 1.0/exposure );
    hLightCurve->Fill( holdTime );
  }
  std::cout << "Randomize SRC (time, energy) ---- done" << std::endl;

  /*
  // randamize SRC (time, energy) -----------------------------------
  for(int i=0; i<countsSrc; i++){
    double srcRandTime = gRandom->Rndm() * v_inteLc[ lcBins-1 ];
    double srcRandInte = gRandom->Rndm() * v_inteSrc[ specBins-1 ];
    double holdTime, holdEnergy;

    for(int j=0; j<lcBins; j++){
      double Y1 = v_inteLc[ j ];
      double Y2 = v_inteLc[ j+1 ];
      if(Y1 < srcRandTime && srcRandTime <= Y2){
	double X1 = v_inteStepLc[ j ];
	double X2 = v_inteStepLc[ j+1 ];
	double A = linearConst(X1, Y1, X2, Y2);
	double B = linearIndex(X1, Y1, X2, Y2);
	double v = linearInvValue(srcRandTime, A, B);
	holdTime = v;
	break;
      }
    }


    for(int j=0; j<specBins; j++){
      double Y1 = v_inteSrc[ j ];
      double Y2 = v_inteSrc[ j+1 ];
      if(Y1 < srcRandInte && srcRandInte <= Y2){
	double X1 = v_inteStepSrc[ j ];
	double X2 = v_inteStepSrc[ j+1 ];
	//double A = linearConst(X1, Y1, X2, Y2);
	//double B = linearIndex(X1, Y1, X2, Y2);
	//double v = linearInvValue(srcRandInte, A, B);
	double A = powConst(X1, Y1, X2, Y2);
	double B = powIndex(X1, Y1, X2, Y2);
	double v = powInvValue(srcRandInte, A, B);
	if( !isnan(v) || !isinf(v) || v<0.1 ){
	  A = linearConst(X1, Y1, X2, Y2);
	  B = linearIndex(X1, Y1, X2, Y2);
	  v = linearInvValue(srcRandInte, A, B);
	}
	holdEnergy = v;
	break;
      }
    }

    v_TandE.push_back(std::pair<double, double>(holdTime, holdEnergy));
    hCountSrc->Fill( holdEnergy, 1.0/exposure );
    hCount->Fill( holdEnergy, 1.0/exposure );
    hCountSrcQdp->Fill( holdEnergy, 1.0/exposure );
    hCountQdp->Fill( holdEnergy, 1.0/exposure );
    hLightCurve->Fill( holdTime );
  }
  */
  
  // sort by time --------------------
  std::sort(v_TandE.begin(), v_TandE.end());
  std::cout << "Event Number(raw): " << v_TandE.size() << std::endl;

  //--------------------------------------------------------------------
  // pile-up
  double pileConstA = 2.0*3.0*3.0*1.0e-14;
  double pileConstB = 2.0*6.5*6.5*1.0e-14;
  double pileConstC = pow(10.0,7);
  const int forNumber = counts + overflowCounts;
  
  for(int i=0; i<forNumber; i++){
    double mainPeakTime = v_TandE[i].first; // [sec]
    double mainPeakEnergy = v_TandE[i].second; // [keV]
    
    if(i%10000==0) cout << "(" << i <<"/" << forNumber << ")" << endl;
    
    if( i != 0 ){ // SAKANABORU
      int j = i - 1;
      double anotherPeakTime = v_TandE[ j ].first;
      double diffTime = mainPeakTime - anotherPeakTime;
      while( diffTime < 0.0000106 ){
	if( diffTime < 0.000002 ){ // zone 2
	  double anotherPeakEnergy = v_TandE[ j ].second;
	  mainPeakEnergy += (anotherPeakEnergy*(4.5/4.0*exp(-diffTime*diffTime/pileConstB) - 0.125));
	}
	else{ // zoon 3
	  double anotherPeakEnergy = v_TandE[ j ].second;
	  mainPeakEnergy += (anotherPeakEnergy*(0.125/80.0*(diffTime*pileConstC-25.0)-0.125));
	}
	j--;
	if( j < 0 ){ break; }
	anotherPeakTime = v_TandE[ j ].first;
	diffTime = mainPeakTime - anotherPeakTime;
      }
    }

    if( i != forNumber - 1 ){ // SUSUMU
      int j = i + 1;
      double anotherPeakTime = v_TandE[ j ].first;
      double diffTime = anotherPeakTime - mainPeakTime;
      while( diffTime < 0.000001 ){ // zone 1
        double anotherPeakEnergy = v_TandE[j].second;
        mainPeakEnergy += (anotherPeakEnergy*exp(-diffTime*diffTime/pileConstA));
	j++;
	if( j > forNumber - 1 ){ break; }
	anotherPeakTime = v_TandE[ j ].first;
	diffTime = anotherPeakTime - mainPeakTime;
      }
    }
    v_TandPileE.push_back(std::pair<double, double>(mainPeakTime, mainPeakEnergy));
  }
  std::cout << "Pileup calculation is done." << std::endl;

  //------------------------------------------------------------------
  // dead time -------------------------------------------------------
  double time = 0.0;
  double trigTiming = 0.000001; // 1 [us]
  const double LD2 = 0.0000024; //  2.4 [us]
  const double overUD = 0.0000256; // 25.6 [us]
  const double normalLD = 0.0000126; // 12.6 [us]

  for(int i=0; i<forNumber; i++){
    if(i%10000==0) cout << "(" << i <<"/" << forNumber << ")" << endl;

    double mainPeakTime = v_TandPileE[i].first; // [sec]
    double trigPeakTime = mainPeakTime - trigTiming; // [sec]
    double mainPeakEnergy = v_TandPileE[i].second; // [keV]

    if( time <= mainPeakTime && trigPeakTime < time && mainPeakEnergy > threshEnergy){ // deadtime expand
      countDeadTime += mainPeakTime - time + (LD2*0.5);
      time = mainPeakTime + (LD2*0.5);
      continue;
    }
    
    if( mainPeakEnergy > threshEneUpp ){ // UD event
      if( trigPeakTime < time ){
	countDeadTime += (trigPeakTime + overUD) - time;
      }
      else{
	countDeadTime += overUD;
      }
      time = trigPeakTime + overUD;
      continue;
    }

    if( mainPeakEnergy < threshEnergy ){ // not over LD
      continue;
    }
    
    if( trigPeakTime < time ) continue;  // the event in deadtime
    
    
    // search maximum (peak) event ------------------
    std::vector<double> peakEnergy, peakTime;
    double trigRegion = trigPeakTime + LD2;
    double peakTimeInThisRegion = trigPeakTime;

    int j = i + 1;
    int jNum = 0;
    if( j < forNumber ){
      while( v_TandPileE[j].first < trigRegion ){
	peakTime.push_back( v_TandPileE[j].first );
	peakEnergy.push_back( v_TandPileE[j].second );
	jNum++;
	j++;
	i++;
	if( j >= forNumber )break;
      }
    }

    int breakFlag = 0;
    for(int k=0; k<jNum; k++){
      if( peakEnergy[ k ] > mainPeakEnergy ){
	mainPeakEnergy = peakEnergy[ k ];
	if( peakTime[ k ] > trigRegion-(trigTiming*0.1) ){ // not Peak during LD2
	  breakFlag = 1;
	  break;
	}
      }
    }
    if( breakFlag == 1 ){ // not PEAK during LD2
      countDeadTime += normalLD;
      time = trigPeakTime + normalLD; // normal one event : 12.6 um
      continue;
    }

    if( mainPeakEnergy > threshEneUpp){ // over UD
      countDeadTime += overUD;
      time = trigPeakTime + overUD;
      continue;
    }

    // fill energy ---------------------------
    hCountDead->Fill( mainPeakEnergy );
    hCountDeadQdp->Fill( mainPeakEnergy );
    countDeadTime += normalLD;
    time = trigPeakTime + normalLD; // normal one event : 12.6 um
    countLiveNumber++;
  }
  std::cout << "Deadtime calculation is done." << std::endl;
  std::cout << "Event Number(pileuped): " << countLiveNumber << std::endl;
  std::cout << "Exposure(sec): " << exposure << std::endl;
  std::cout << "Dead time(sec): " << countDeadTime << std::endl;
}


void simPileup::runNoPileup(){
  if( lcBins == 0 || specBins == 0 ){
    std::cout << "ERROR: LC and/or SPEC have not been set." << std::endl;
    exit(1);
  }
  
  // count rate ------------------------------
  countsBgd = (int)(exposure * countsBgd);
  countsSrc = counts - countsBgd;
  if( counts <= countsBgd ){
    std::cout << "ERROR: Incident counts is too less than BGD counts." << std::endl;
    exit(1);
  }
  
  // randamize BGD (time, energy) -----------------------------------
  for(int i=0; i<countsBgd; i++){
    double bgdTime = gRandom->Rndm() * exposure;
    double bgdRandInte = gRandom->Rndm() * v_inteBgd[ specBins-1 ];
    
    for(int j=0; j<specBins; j++){
      double Y1 = v_inteBgd[ j ];
      double Y2 = v_inteBgd[ j+1 ];
      if(Y1 < bgdRandInte && bgdRandInte <= Y2){
	double X1 = v_inteStepBgd[ j ];
	double X2 = v_inteStepBgd[ j+1 ];
	//double A = linearConst(X1, Y1, X2, Y2);
	//double B = linearIndex(X1, Y1, X2, Y2);
	//double v = linearInvValue(bgdRandInte, A, B);
	double A = powConst(X1, Y1, X2, Y2);
	double B = powIndex(X1, Y1, X2, Y2);
	double v = powInvValue(bgdRandInte, A, B);
	if( !isnan(v) || !isinf(v) || v<0.1 ){
	  A = linearConst(X1, Y1, X2, Y2);
	  B = linearIndex(X1, Y1, X2, Y2);
	  v = linearInvValue(bgdRandInte, A, B);// [keV]
	}
	v_TandE.push_back(std::pair<double, double>(bgdTime, v));
	hCountBgd->Fill( v, 1.0/exposure );
	hCount->Fill( v, 1.0/exposure );
	hCountBgdQdp->Fill( v, 1.0/exposure );
	hCountQdp->Fill( v, 1.0/exposure );
	hLightCurve->Fill( bgdTime );
	break;
      }
    }
  }

  // randamize SRC (time, energy) -----------------------------------
  for(int i=0; i<countsSrc; i++){
    double srcRandTime = gRandom->Rndm() * v_inteLc[ lcBins-1 ];
    double srcRandInte = gRandom->Rndm() * v_inteSrc[ specBins-1 ];
    double holdTime, holdEnergy;

    for(int j=0; j<lcBins; j++){
      double Y1 = v_inteLc[ j ];
      double Y2 = v_inteLc[ j+1 ];
      if(Y1 < srcRandTime && srcRandTime <= Y2){
	double X1 = v_inteStepLc[ j ];
	double X2 = v_inteStepLc[ j+1 ];
	double A = linearConst(X1, Y1, X2, Y2);
	double B = linearIndex(X1, Y1, X2, Y2);
	double v = linearInvValue(srcRandTime, A, B);
	holdTime = v;
	break;
      }
    }


    for(int j=0; j<specBins; j++){
      double Y1 = v_inteSrc[ j ];
      double Y2 = v_inteSrc[ j+1 ];
      if(Y1 < srcRandInte && srcRandInte <= Y2){
	double X1 = v_inteStepSrc[ j ];
	double X2 = v_inteStepSrc[ j+1 ];
	//double A = linearConst(X1, Y1, X2, Y2);
	//double B = linearIndex(X1, Y1, X2, Y2);
	//double v = linearInvValue(srcRandInte, A, B);
	double A = powConst(X1, Y1, X2, Y2);
	double B = powIndex(X1, Y1, X2, Y2);
	double v = powInvValue(srcRandInte, A, B);
	if( !isnan(v) || !isinf(v) || v<0.1 ){
	  A = linearConst(X1, Y1, X2, Y2);
	  B = linearIndex(X1, Y1, X2, Y2);
	  v = linearInvValue(srcRandInte, A, B);
	}
	holdEnergy = v;
	break;
      }
    }

    v_TandE.push_back(std::pair<double, double>(holdTime, holdEnergy));
    hCountSrc->Fill( holdEnergy, 1.0/exposure );
    hCount->Fill( holdEnergy, 1.0/exposure );
    hCountSrcQdp->Fill( holdEnergy, 1.0/exposure );
    hCountQdp->Fill( holdEnergy, 1.0/exposure );
    hLightCurve->Fill( holdTime );
  }

  
  // sort by time --------------------
  std::sort(v_TandE.begin(), v_TandE.end());
  std::cout << "Event Number(raw): " << v_TandE.size() << std::endl;
  /*
  //--------------------------------------------------------------------
  // pile-up
  double pileConstA = 2.0*3.0*3.0*1.0e-14;
  double pileConstB = 2.0*6.5*6.5*1.0e-14;
  double pileConstC = pow(10.0,7);
  
  for(int i=0; i<counts + overflowCounts; i=i+1){
    double mainPeakTime = v_TandE[i].first; // [sec]
    double mainPeakEnergy = v_TandE[i].second; // [keV]
    
    if(i%10000==0) cout << "(" << i <<"/" << counts + overflowCounts << ")" << endl;
    
    int start_j = i-nearEventNum;
    int end_j = i+nearEventNum;
    if(start_j < 0) start_j = 0;
    if(end_j > counts + overflowCounts) end_j = counts + overflowCounts;
    
   for(int j=start_j; j<end_j; j=j+1){
      if(i==j)continue;
      
      double anotherPeakTime = v_TandE[j].first;
      if(mainPeakTime <= anotherPeakTime && anotherPeakTime < mainPeakTime+0.000001){ // zone 1
	double anotherPeakEnergy = v_TandE[j].second;
	mainPeakEnergy += (anotherPeakEnergy*exp(-(anotherPeakTime-mainPeakTime)*(anotherPeakTime-mainPeakTime)/pileConstA));
      }
      else if(mainPeakTime-0.000002 <= anotherPeakTime && anotherPeakTime < mainPeakTime){ // zone 2
	double anotherPeakEnergy = v_TandE[j].second;
	mainPeakEnergy += (anotherPeakEnergy*(4.5/4.0*exp(-(anotherPeakTime-mainPeakTime)*(anotherPeakTime-mainPeakTime)/pileConstB) - 0.125));
      }
      else if(mainPeakTime-0.0000106 <= anotherPeakTime && anotherPeakTime < mainPeakTime-0.000002){ // zone 3
	double anotherPeakEnergy = v_TandE[j].second;
	mainPeakEnergy += (anotherPeakEnergy*(0.125/80.0*((mainPeakTime-anotherPeakTime)*pileConstC-25.0)-0.125));
      }
   }   
   v_TandPileE.push_back(std::pair<double, double>(mainPeakTime, mainPeakEnergy));
  }
  std::cout << "Pileup calculation is done." << std::endl;

  //------------------------------------------------------------------
  // dead time -------------------------------------------------------
  double time = 0.0;
  double trigTiming = 0.000001; // 1 [us]
  const int forNumber = counts + overflowCounts;
  const double LD2 = 0.0000024; //  2.4 [us]
  const double overUD = 0.0000256; // 25.6 [us]
  const double normalLD = 0.0000126; // 12.6 [us]

  for(int i=0; i<forNumber; i=i+1){
    if(i%10000==0) cout << "(" << i <<"/" << forNumber << ")" << endl;

    double mainPeakTime = v_TandPileE[i].first; // [sec]
    double trigPeakTime = mainPeakTime - trigTiming; // [sec]
    double mainPeakEnergy = v_TandPileE[i].second; // [keV]

    //if( mainPeakTime > time-trigTiming && mainPeakTime < time && mainPeakEnergy > threshEnergy){ // when any triggers are active.
    //  time += mainPeakTime - time;
    //  countDeadTime += mainPeakTime - time;
    //  continue;
    //}  
    if( time <= mainPeakTime && mainPeakTime < time+trigTiming && mainPeakEnergy > threshEnergy){ // when any triggers are active.
      time += mainPeakTime - time;
      countDeadTime += mainPeakTime - time;
      continue;
    }
    if( trigPeakTime < time && mainPeakEnergy > threshEneUpp){ // event is in dead time      
      //time += (trigPeakTime + overUD) - time;
      //countDeadTime += (trigPeakTime + overUD) - time;
      time += (trigPeakTime + normalLD) - time;
      countDeadTime += (trigPeakTime + normalLD) - time;
      continue;
    }
    else if( trigPeakTime < time ) continue;            // event is in dead time
    if( mainPeakEnergy < threshEnergy ){           // not over LD
      time += mainPeakTime-time;
      continue;
    }
    if( mainPeakEnergy > threshEneUpp){            // over UD
      time += (trigPeakTime-time) + overUD;
      countDeadTime += overUD;
      continue;
    }
    

    // search maximum (peak) event ------------------
    std::vector<double> peakEnergy, peakTime;
    double trigRegion = trigPeakTime + LD2;

    int j = i + 1;
    int jNum = 0;
    if( j < forNumber ){
      while( v_TandPileE[j].first < trigRegion ){
	peakTime.push_back( v_TandPileE[j].first );
	peakEnergy.push_back( v_TandPileE[j].second );
	jNum++;
	j++;
      }
    }

    int breakFlag = 0;
    for(int k=0; k<jNum; k=k+1){
      if( peakEnergy[ k ] > mainPeakEnergy ) mainPeakEnergy = peakEnergy[ k ];
      if( peakTime[ k ] > trigRegion-(trigTiming*0.1) ){ // not Peak during LD2
	breakFlag = 1;
	break;
      }
    }
    if( breakFlag == 1 ){ // not PEAK during LD2
      time += (trigPeakTime-time) + normalLD; // normal one event : 12.6 um
      countDeadTime += normalLD;
      continue;
    }
    
    if( mainPeakEnergy > threshEneUpp){            // over UD
      time += (trigPeakTime-time) + overUD;
      countDeadTime += overUD;
      continue;
    }

    // fill energy ---------------------------
    hCountDead->Fill( mainPeakEnergy );
    hCountDeadQdp->Fill( mainPeakEnergy );
    time += (trigPeakTime-time) + normalLD; // normal one event : 12.6 um
    countDeadTime += normalLD;
    countLiveNumber++;
  }
  std::cout << "Deadtime calculation is done." << std::endl;
  std::cout << "Event Number(pileuped): " << countLiveNumber << std::endl;
  std::cout << "Exposure(sec): " << exposure << std::endl;
  std::cout << "Dead time(sec): " << countDeadTime << std::endl;
*/
}


double simPileup::getDeadtime(){
  return countDeadTime;
}

double simPileup::getLiveCount(){
  return countLiveNumber;
}


sli::mdarray_double simPileup::getOriginMatrix(){
  sli::mdarray_double matrix;
  matrix.init( false );
  matrix.resize_2d(hCount->GetXaxis()->GetNbins(), 1);
  for(int i=1; i<=hCount->GetXaxis()->GetNbins(); i++){
    double val = hCount->GetBinContent( i );
    matrix(i-1, 0) = val;
  }
  return matrix;
}


sli::mdarray_double simPileup::getSrcMatrix(){
  sli::mdarray_double matrix;
  matrix.init( false );
  matrix.resize_2d(hCountSrc->GetXaxis()->GetNbins(), 1);
  for(int i=1; i<=hCountSrc->GetXaxis()->GetNbins(); i++){
    double val = hCountSrc->GetBinContent( i );
    matrix(i-1, 0) = val;
  }
  return matrix;
}


sli::mdarray_double simPileup::getBgdMatrix(){
  sli::mdarray_double matrix;
  matrix.init( false );
  matrix.resize_2d(hCountBgd->GetXaxis()->GetNbins(), 1);
  for(int i=1; i<=hCountBgd->GetXaxis()->GetNbins(); i++){
    double val = hCountBgd->GetBinContent( i );
    matrix(i-1, 0) = val;
  }
  return matrix;
}


sli::mdarray_double simPileup::getDeadMatrix(){
  sli::mdarray_double matrix;
  matrix.init( false );
  matrix.resize_2d(hCountDead->GetXaxis()->GetNbins(), 1);
  for(int i=1; i<=hCountDead->GetXaxis()->GetNbins(); i++){
    double val = hCountDead->GetBinContent( i );
    matrix(i-1, 0) = val;
  }
  return matrix;
}

void simPileup::outputLightCurve(){
  ofstream output_lcQdp;
  output_lcQdp.open("output_lightcurve.qdp");
  output_lcQdp << "read serr 1 2" << std::endl;
  output_lcQdp << "line on" << std::endl;
  output_lcQdp << "lw 5" << std::endl;
  output_lcQdp << "la x Time (sec)" << std::endl;
  for(int i=1; i<hLightCurve->GetXaxis()->GetNbins()+1; i++){
    double binCentor = hLightCurve->GetXaxis()->GetBinCenter( i );
    double binWidth = hLightCurve->GetXaxis()->GetBinWidth( i )*0.5;
    double val = hLightCurve->GetBinContent(i);
    output_lcQdp << binCentor <<" "<< binWidth <<" "<< val <<" "<< sqrt( val ) <<std::endl;
  }
  output_lcQdp.close(); 
}


void simPileup::outputOriginSpec(){
  ofstream output_Qdp;
  output_Qdp.open("output_originSpec.qdp");
  output_Qdp << "read serr 1 2" << std::endl;
  output_Qdp << "line off" << std::endl;
  output_Qdp << "log" << std::endl;
  output_Qdp << "lw 5" << std::endl;
  output_Qdp << "la x Energy (keV)" << std::endl;
  output_Qdp << "la y normalized count/s/keV" << std::endl;
  for(int i=1; i<hCountQdp->GetXaxis()->GetNbins()+1; i++){
    double binCentor = hCountQdp->GetXaxis()->GetBinCenter( i );
    double binWidth = hCountQdp->GetXaxis()->GetBinWidth( i )*0.5;
    double val = hCountQdp->GetBinContent(i);
    output_Qdp << binCentor <<" "<< binWidth <<" "<< val/(binWidth*2.0) <<" "<< sqrt( val )/(binWidth*2.0) <<std::endl;
  }
  output_Qdp.close(); 
}

void simPileup::outputSrcSpec(){
  ofstream output_Qdp;
  output_Qdp.open("output_srcSpec.qdp");
  output_Qdp << "read serr 1 2" << std::endl;
  output_Qdp << "line off" << std::endl;
  output_Qdp << "log" << std::endl;
  output_Qdp << "lw 5" << std::endl;
  output_Qdp << "la x Energy (keV)" << std::endl;
  output_Qdp << "la y normalized count/s/keV" << std::endl;
  for(int i=1; i<hCountSrcQdp->GetXaxis()->GetNbins()+1; i++){
    double binCentor = hCountSrcQdp->GetXaxis()->GetBinCenter( i );
    double binWidth = hCountSrcQdp->GetXaxis()->GetBinWidth( i )*0.5;
    double val = hCountSrcQdp->GetBinContent(i);
    output_Qdp << binCentor <<" "<< binWidth <<" "<< val/(binWidth*2.0) <<" "<< sqrt( val )/(binWidth*2.0) <<std::endl;
  }
  output_Qdp.close(); 
}

void simPileup::outputBgdSpec(){
  ofstream output_Qdp;
  output_Qdp.open("output_bgdSpec.qdp");
  output_Qdp << "read serr 1 2" << std::endl;
  output_Qdp << "line off" << std::endl;
  output_Qdp << "log" << std::endl;
  output_Qdp << "lw 5" << std::endl;
  output_Qdp << "la x Energy (keV)" << std::endl;
  output_Qdp << "la y normalized count/s/keV" << std::endl;
  for(int i=1; i<hCountBgdQdp->GetXaxis()->GetNbins()+1; i++){
    double binCentor = hCountBgdQdp->GetXaxis()->GetBinCenter( i );
    double binWidth = hCountBgdQdp->GetXaxis()->GetBinWidth( i )*0.5;
    double val = hCountBgdQdp->GetBinContent(i);
    output_Qdp << binCentor <<" "<< binWidth <<" "<< val/(binWidth*2.0) <<" "<< sqrt( val )/(binWidth*2.0) <<std::endl;
  }
  output_Qdp.close(); 
}

void simPileup::outputDeadSpec(){
  ofstream output_Qdp;
  output_Qdp.open("output_deadSpec.qdp");
  output_Qdp << "read serr 1 2" << std::endl;
  output_Qdp << "line off" << std::endl;
  output_Qdp << "log" << std::endl;
  output_Qdp << "lw 5" << std::endl;
  output_Qdp << "la x Energy (keV)" << std::endl;
  output_Qdp << "la y normalized count/s/keV" << std::endl;
  for(int i=1; i<hCountDeadQdp->GetXaxis()->GetNbins()+1; i++){
    double binCentor = hCountDeadQdp->GetXaxis()->GetBinCenter( i );
    double binWidth = hCountDeadQdp->GetXaxis()->GetBinWidth( i )*0.5;
    double val = hCountDeadQdp->GetBinContent(i);
    output_Qdp << binCentor <<" "<< binWidth <<" "<< val/(binWidth*2.0)/(1.0-countDeadTime) <<" "<< sqrt( val )/(binWidth*2.0)/(1.0-countDeadTime) <<std::endl;
  }
  output_Qdp.close(); 
}
