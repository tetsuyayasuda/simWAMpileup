//--------------------------------------------
// readText.hh
//
// 
//--------------------------------------------
// 2012-08-11 yasuda ver0.0

// Usage:
//  #include "readText.hh"
//  readText* test = new readText("inputFile.dat");
//  test->setColumn3();  
//           // when inputFile.dat has three columns. (ex: light curve)
//  test->setColumn4();  
//           //when inputFile.dat has four columns. (ex: spectrum)
//  cout << test->getColumn2(10) << endl;  
//           // if you want to read the value of second column of 10th line.
//  cout << test->getColumn3(10, 2) << endl;  
//           // if you want to read the value of second column of 10th line.
//

#ifndef READTEXT_HH
#define READTEXT_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

class readText
{

public:
  readText();
  void setFile(char* input);
  virtual void setColumn1();
  virtual void setColumn2();
  virtual void setColumn3();
  virtual void setColumn4();
  
  float getColumn1(int eventNum, int columnNum=1);
  float getColumn2(int eventNum, int columnNum=2);
  float getColumn3(int eventNum, int columnNum=3);
  float getColumn4(int eventNum, int columnNum=4);
 
  int getColumnNum();
  int getEventNum();

  //private:
  char* fileName;
  ifstream inputFile;
  int columnNumber;
  int eventNumber;

  vector<float> column1;
  vector<float> column2;
  vector<float> column3;
  vector<float> column4;

};

#endif
