//--------------------------------------------
//     simPileup.hh
//--------------------------------------------
// 2013-09-13 yasuda

#ifndef SIMPILEUP_HH
#define SIMPILEUP_HH
#include <string>
#include <vector>

#include <sli/stdstreamio.h>
#include <sli/tstring.h>
#include <sli/fitscc.h>
#include <sli/mdarray.h>
#include <sli/mdarray_statistics.h>
#include <sli/mdarray_math.h>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TLinearMinimizer.h>
#include <TH1D.h>

class simPileup
{
public:
  simPileup();
  ~simPileup();
  void setEnergyRange(double eBounds);
  void setCountRate(int countRate);
  void preparation();
  void inteLightcurve();
  void inteLightcurve(char* lcName);
  void inteSrc(sli::mdarray_double matrix);
  void inteBgd(sli::mdarray_double matrix);
  void setOverflowRate(int overflowRate, double liveTime);
  void runPileup();
  void runNoPileup();
  double getDeadtime();
  double getLiveCount();
  
  sli::mdarray_double getOriginMatrix();
  sli::mdarray_double getSrcMatrix();
  sli::mdarray_double getBgdMatrix();
  sli::mdarray_double getDeadMatrix();

  void outputLightCurve();
  void outputOriginSpec();
  void outputSrcSpec();
  void outputBgdSpec();
  void outputDeadSpec();

private:
  int nearEventNum;
  double threshEnergy; // [keV]
  double threshEneUpp; // [keV]
  std::vector<double> v_energyBounds, v_energyBoundsQdp;
  int numBounds;
  int counts, countsBgd, countsSrc;
  std::vector<double> v_inteLc, v_inteStepLc;
  std::vector<double> v_inteSrc, v_inteStepSrc;
  std::vector<double> v_inteBgd, v_inteStepBgd;
  std::vector<std::pair<double, double> > v_TandE;
  std::vector<std::pair<double, double> > v_TandPileE;
  std::vector<std::pair<double, double> >::iterator it;
  int lcBins, specBins;
  double exposure;
  double countDeadTime;
  int countLiveNumber;
  int separateStepLc, separateStepSpec;
  int overflowCounts;

  TH1D *hCount, *hCountQdp; // raw events (SRC+BGD)
  TH1D *hCountSrc, *hCountSrcQdp; // raw events (SRC)
  TH1D *hCountBgd, *hCountBgdQdp; // raw events (BGD)
  TH1D *hCountDead, *hCountDeadQdp; // pileup events (pileupSRC + pileupBGD)
  TH1D *hLightCurve;

  sli::mdarray_double matrixSrc;
  sli::mdarray_double matrixBgd;


};

#endif
