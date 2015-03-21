//--------------------------------------------
// readText.cc
//
// 
//--------------------------------------------
// 2012-08-11 yasuda Ver0.0


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "readText.hh"

using namespace std;

readText::readText()
{
}

void readText::setFile(char* input)
{
  // initialize --------------------------
  fileName = input;  
  columnNumber = 0;
  eventNumber = 0;
}

int readText::getColumnNum(){
  return columnNumber;
}

int readText::getEventNum(){
  return eventNumber;
}

void readText::setColumn1()
{
  // file open ---------------------------
  inputFile.open(fileName);
  if( ! inputFile.good() ){
    cout << "--------------  ERROR  ---------------" << endl;
    cout << "Could not open: " << fileName << " !!!" << endl;
    exit(1);
  }

  // initialize --------------------------
  columnNumber = 1;

  // take Value --------------------------
  while( ! inputFile.eof() )
    {
      if(column1.size() <= eventNumber ){
	column1.resize( eventNumber + 1);
      }
      inputFile >> column1[eventNumber];
      eventNumber = eventNumber + 1;
    }
  eventNumber = eventNumber -1;
  
  // comment -----------------------------
  cout << "Column Number = " << columnNumber << endl;
  cout << "Event Number = " << eventNumber << endl;
  
  // close file --------------------------
  inputFile.close();
}


void readText::setColumn2()
{
  // file open ---------------------------
  inputFile.open(fileName);
  if( ! inputFile.good() ){
    cout << "--------------  ERROR  ---------------" << endl;
    cout << "Could not open: " << fileName << " !!!" << endl;
    exit(1);
  }

  // initialize --------------------------
  columnNumber = 2;

  // take Value --------------------------
  while( ! inputFile.eof() )
    {
      if(column1.size() <= eventNumber ){
	column1.resize( eventNumber + 1);
	column2.resize( eventNumber + 1);
      }
      inputFile >> column1[eventNumber] >> column2[eventNumber];
      eventNumber = eventNumber + 1;
    }
  eventNumber = eventNumber -1;
  
  // comment -----------------------------
  cout << "Column Number = " << columnNumber << endl;
  cout << "Event Number = " << eventNumber << endl;

  // close file --------------------------
  inputFile.close();
}


void readText::setColumn3()
{
  // file open ---------------------------
  inputFile.open(fileName);
  if( ! inputFile.good() ){
    cout << "--------------  ERROR  ---------------" << endl;
    cout << "Could not open: " << fileName << " !!!" << endl;
    exit(1);
  }

  // initialize --------------------------
  columnNumber = 3;

  // take Value --------------------------
  while( ! inputFile.eof() )
    {
      if(column1.size() <= eventNumber ){
	column1.resize( eventNumber + 1);
	column2.resize( eventNumber + 1);
	column3.resize( eventNumber + 1);
      }
      inputFile >> column1[eventNumber] >> column2[eventNumber] >> column3[eventNumber];
      eventNumber = eventNumber + 1;
    }
  eventNumber = eventNumber -1;
  
  // comment -----------------------------
  cout << "Column Number = " << columnNumber << endl;
  cout << "Event Number = " << eventNumber << endl;

  // close file --------------------------
  inputFile.close();
}

 
void readText::setColumn4()
{
  // file open ---------------------------
  inputFile.open(fileName);
  if( ! inputFile.good() ){
    cout << "--------------  ERROR  ---------------" << endl;
    cout << "Could not open: " << fileName << " !!!" << endl;
    exit(1);
  }

  // initialize --------------------------
  columnNumber = 4;

  // take Value --------------------------
  while( ! inputFile.eof() )
    {
      if(column1.size() <= eventNumber ){
	column1.resize( eventNumber + 1);
	column2.resize( eventNumber + 1);
	column3.resize( eventNumber + 1);
	column4.resize( eventNumber + 1);
      }
      inputFile >> column1[eventNumber] >> column2[eventNumber] >> column3[eventNumber] >> column4[eventNumber];
      eventNumber = eventNumber + 1;
    }
  eventNumber = eventNumber -1;
  
  // comment -----------------------------
  cout << "Column Number = " << columnNumber << endl;
  cout << "Event Number = " << eventNumber << endl;

  // close file --------------------------
  inputFile.close();
}


float readText::getColumn1(int eNum=0, int cNum){
  // define return value ------------------
  float value;
  
  // find invalid value -------------------
  if(eNum < 0 || eNum > eventNumber || cNum < 1 || cNum > columnNumber){
    cout << "----------- ERROR -----------" << endl;
    cout << "invalid parametars [getColumn1]" << endl;
    exit(1);
  }
  
  // get value ----------------------------
  value = column1[eNum];

  return value;
}


float readText::getColumn2(int eNum=0, int cNum){
  // define return value ------------------
  float value;

  // find invalid value -------------------
  if(eNum < 0 || eNum > eventNumber || cNum < 1 || cNum > columnNumber){
    cout << "----------- ERROR -----------" << endl;
    cout << "invalid parametars [getColumn2]" << endl;
    exit(1);
  }
  
  // get value ----------------------------
  if(cNum == 1) value = column1[eNum];
  else value = column2[eNum];

  return value;
}


float readText::getColumn3(int eNum=0, int cNum){
  // define return value ------------------
  float value;

  // find invalid value -------------------
  if(eNum < 0 || eNum > eventNumber || cNum < 1 || cNum > columnNumber){
    cout << "----------- ERROR -----------" << endl;
    cout << "invalid parametars [getColumn3]" << endl;
    exit(1);
  }
  
  // get value ----------------------------
  if(cNum == 1) value = column1[eNum];
  else if(cNum == 2) value = column2[eNum];
  else value = column3[eNum];

  return value;
}


float readText::getColumn4(int eNum=0, int cNum){
  // define return value ------------------
  float value;

  // find invalid value -------------------
  if(eNum < 0 || eNum > eventNumber || cNum < 1 || cNum > columnNumber){
    cout << "----------- ERROR -----------" << endl;
    cout << "invalid parametars [getColumn4]" << endl;
    exit(1);
  }
  
  // get value ----------------------------
  if(cNum == 1) value = column1[eNum];
  else if(cNum == 2) value = column2[eNum];
  else if(cNum == 3) value = column3[eNum];
  else value = column4[eNum];

  return value;
}


