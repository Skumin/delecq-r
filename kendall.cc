#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>

#ifdef SOLARIS_OLD      /* 'if' might be overkill, but just to be minimal... */
/*                               10 May 2010 [rickr] */
#include "machdep.h"
#else
#include <stdint.h>
#endif
  
  
  /*-------------------------------------------------------------------------*/
  /* Sorts in place, returns the bubble sort distance between the input array
* and the sorted array.
*/
  
  static int insertionSort(double *arr, int len)
{
  int maxJ, i,j , swapCount = 0;
  
  /* printf("enter insertionSort len=%d\n",len) ; */
    
    if(len < 2) { return 0; }
  
  maxJ = len - 1;
  for(i = len - 2; i >= 0; --i) {
    double  val = arr[i];
    for(j=i; j < maxJ && arr[j + 1] < val; ++j) {
      arr[j] = arr[j + 1];
    }
    
    arr[j] = val;
    swapCount += (j - i);
  }
  
  return swapCount;
  }

/*-------------------------------------------------------------------------*/
  
  static int merge(double *from, double *to, int middle, int len)
{
  int bufIndex, leftLen, rightLen , swaps ;
  double *left , *right;
  
  /* printf("enter merge\n") ; */
    
    bufIndex = 0;
  swaps = 0;
  
  left = from;
  right = from + middle;
  rightLen = len - middle;
  leftLen = middle;
  
  while(leftLen && rightLen) {
    if(right[0] < left[0]) {
      to[bufIndex] = right[0];
      swaps += leftLen;
      rightLen--;
      right++;
    } else {
      to[bufIndex] = left[0];
      leftLen--;
      left++;
    }
    bufIndex++;
  }
  
  if(leftLen) {
    #pragma omp critical (MEMCPY)
    memcpy(to + bufIndex, left, leftLen * sizeof(double));
  } else if(rightLen) {
    #pragma omp critical (MEMCPY)
    memcpy(to + bufIndex, right, rightLen * sizeof(double));
  }
  
  return swaps;
  }

/*-------------------------------------------------------------------------*/
  /* Sorts in place, returns the bubble sort distance between the input array
* and the sorted array.
*/
  
  static int mergeSort(double *x, double *buf, int len)
{
  int swaps , half ;
  
  /* printf("enter mergeSort\n") ; */
    
    if(len < 10) {
      return insertionSort(x, len);
    }
  
  swaps = 0;
  
  if(len < 2) { return 0; }
  
  half = len / 2;
  swaps += mergeSort(x, buf, half);
  swaps += mergeSort(x + half, buf + half, len - half);
  swaps += merge(x, buf, half, len);
  
  #pragma omp critical (MEMCPY)
  memcpy(x, buf, len * sizeof(double));
  return swaps;
  }

/*-------------------------------------------------------------------------*/
  
  static int getMs(double *data, int len)  /* Assumes data is sorted */
  {
    int Ms = 0, tieCount = 0 , i ;
    
    /* printf("enter getMs\n") ; */
      
      for(i = 1; i < len; i++) {
        if(data[i] == data[i-1]) {
          tieCount++;
        } else if(tieCount) {
          Ms += (tieCount * (tieCount + 1)) / 2;
          tieCount = 0;
        }
      }
    if(tieCount) {
      Ms += (tieCount * (tieCount + 1)) / 2;
    }
    return Ms;
  }

/*-------------------------------------------------------------------------*/
  /* This function calculates the Kendall correlation tau_b.
* The arrays arr1 should be sorted before this call, and arr2 should be
* re-ordered in lockstep.  This can be done by calling
*   qsort_floatfloat(len,arr1,arr2)
* for example.
* Note also that arr1 and arr2 will be modified, so if they need to
* be preserved, do so before calling this function.
*/
  // [[Rcpp::export]]
  double kendallNlogN( double *arr1, double *arr2, int len )
{
  int m1 = 0, m2 = 0, tieCount, swapCount, nPair, s,i ;
  double cor ;
  
  /* printf("enter kendallNlogN\n") ; */
    
    if( len < 2 ) return (double)0 ;
  
  nPair = len * (len - 1) / 2;
  s = nPair;
  
  tieCount = 0;
  for(i = 1; i < len; i++) {
    if(arr1[i - 1] == arr1[i]) {
      tieCount++;
    } else if(tieCount > 0) {
      insertionSort(arr2 + i - tieCount - 1, tieCount + 1);
      m1 += tieCount * (tieCount + 1) / 2;
      s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
      tieCount = 0;
    }
  }
  if(tieCount > 0) {
    insertionSort(arr2 + i - tieCount - 1, tieCount + 1);
    m1 += tieCount * (tieCount + 1) / 2;
    s += getMs(arr2 + i - tieCount - 1, tieCount + 1);
  }
  
  swapCount = mergeSort(arr2, arr1, len);
  
  m2 = getMs(arr2, len);
  s -= (m1 + m2) + 2 * swapCount;
  
  if( m1 < nPair && m2 < nPair )
    cor = s / ( sqrtf((double)(nPair-m1)) * sqrtf((double)(nPair-m2)) ) ;
  else
    cor = 0.0f ;
  
  return cor ;
  }