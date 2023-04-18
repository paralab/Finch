
/**
  @file binUtils.C
  @brief A set of functions for fast binary operations
  @author Rahul S. Sampath, rahul.sampath@gmail.com
  @author Hari Sundar, hsundar@gmail.com
  */

#include <vector>
#include <cassert>
#include "binUtils.h"

namespace binOp {

  unsigned int fastLog2(unsigned int num) {
    if(num) {
      return (binLength(num) - 1);
    } else {
      assert(false);
      return 1; 
    }
  }//end function

  unsigned long long binLength(unsigned long long num) {
    unsigned int len = 1;
    while(num > 1) {
      num = (num >> 1);
      len++;
    }
    return len;
  }//end function

  int toBin(unsigned int num, unsigned int binLen,  std::vector<bool>& numBin) {
    numBin = std::vector<bool>(binLen);
    for(unsigned int i = 0; i < binLen; i++) {
      numBin[i]=0;
    }//end for
    unsigned int pos = binLen -1;
    while(num > 0) {
      numBin[pos] = (num%2);
      num = num/2;  
      pos--;
    }  //end while  
    return 1;
  }//end function

  unsigned int binToDec(unsigned int* numBin, unsigned int binLen) {
    unsigned int res = 0;
    for(unsigned int i = 0; i< binLen; i++) {
      res = (2*res) + numBin[i];
    }
    return res;
  }//end function


  bool isPowerOfTwo(unsigned int n) {
    return (n && (!(n & (n - 1))));
  }

  // compute the next highest power of 2 of 32-bit v
  int getNextHighestPowerOfTwo(unsigned int n) {
    unsigned int v = n;
    assert(v > 0);
    v--;
    v |= (v >> 1);
    v |= (v >> 2);
    v |= (v >> 4);
    v |= (v >> 8);
    v |= (v >> 16);
    v++;
    return v;
  }

  // compute the prev highest power of 2 of 32-bit v
  int getPrevHighestPowerOfTwo(unsigned int n) {
    unsigned int v = n;
    assert(v > 0);
    v--;
    v |= (v >> 1);
    v |= (v >> 2);
    v |= (v >> 4);
    v |= (v >> 8);
    v |= (v >> 16);
    v++;
    return  (v >> 1);
  }

  unsigned int lowestOnePos(unsigned int num)
  {
    assert(num);
    unsigned int pos = 0;
    while (!(num & 1u))
    {
      num >>= 1u;
      pos++;
    }
    return pos;
  }



}//end namespace

