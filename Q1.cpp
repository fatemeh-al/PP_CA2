#include <stdio.h>
#include <x86intrin.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <bits/stdc++.h>

#define ARRAY_SIZE int(pow(2, 20))

using namespace std;

//command:
//g++ Q1.cpp -std=c++11 -msse4.1 && ./a.out

int main(){
    float array[ARRAY_SIZE];
    float a = 100.0, max1 = -1, max2 = -1;
    int maxIndex1 = -1, maxIndex2 = -1;
    timeval startTime, endTime;
    double serialTime, parallelTime;
    float maxValuesSerial[4] = {-1, -1, -1, -1};
    int maxIndexesSerial[4] = {-1, -1, -1, -1};
    __m128 vec;
    __m128 indexes = _mm_set_ps(0.0, 1.0, 2.0, 3.0), increment = _mm_set1_ps(4.0);
    __m128 maxIndexes = _mm_set_ps(-1.0, -1.0, -1.0, -1.0), maxValues = _mm_set_ps(-1.0, -1.0, -1.0, -1.0);
    float tempVec[4], tempIndexVec[4];

    srand((unsigned int)time(NULL));
    for (int i = 0; i < ARRAY_SIZE; i++)
        array[i] = (float(rand())/float((RAND_MAX)) * a);

    //Serial
    gettimeofday(&startTime, NULL);
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        if(array[i] >= maxValuesSerial[0]){
            maxValuesSerial[0] = array[i];
            maxIndexesSerial[0] = i;
        }
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        if(array[i + 1] >= maxValuesSerial[1]){
            maxValuesSerial[1] = array[i + 1];
            maxIndexesSerial[1] = i + 1;
        }
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        if(array[i + 2] >= maxValuesSerial[2]){
            maxValuesSerial[2] = array[i + 2];
            maxIndexesSerial[2] = i + 2;
        }
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        if(array[i + 3] >= maxValuesSerial[3]){
            maxValuesSerial[3] = array[i + 3];
            maxIndexesSerial[3] = i + 3;
        }
    for (int i=0; i < 4; i++) {
        if (maxValuesSerial[i] >= max1) {
            max1 = maxValuesSerial[i];
            maxIndex1 = maxIndexesSerial[i];
        }
    }
    gettimeofday(&endTime, NULL); 
    serialTime = (endTime.tv_sec - startTime.tv_sec) * 1e6; 
    serialTime = (serialTime + (endTime.tv_usec - startTime.tv_usec)) * 1e-6; 

    //Parallel
    gettimeofday(&startTime, NULL);
    for(int i = 0; i < ARRAY_SIZE; i += 4){
        vec = _mm_loadu_ps(&array[i]);
        __m128 gt = _mm_cmpge_ps(vec, maxValues);
        maxIndexes = _mm_blendv_ps(maxIndexes, indexes, gt);
        maxValues = _mm_max_ps(maxValues, vec);

        indexes = _mm_add_ps(indexes, increment);
    }
    _mm_storeu_ps(tempVec, maxValues);
    _mm_storeu_ps(tempIndexVec, maxIndexes);
    for (int i=0; i < 4; i++) {
        if (tempVec[i] >= max2) {
            max2 = tempVec[i];
            maxIndex2 = tempIndexVec[i];
        }
    }
    gettimeofday(&endTime, NULL); 
    parallelTime = (endTime.tv_sec - startTime.tv_sec) * 1e6; 
    parallelTime = (parallelTime + (endTime.tv_usec - startTime.tv_usec)) * 1e-6;

    cout << "Parall Programming - CA2 - Q1 - 810196512 - 810196469" << endl << endl;
    cout << "Serial Time: " << fixed << serialTime << setprecision(6) << " sec" << endl; 
    cout << "Max value: " << max1 << " and its index is: " << maxIndex1 << endl << endl;
    cout << "Parallel Time: " << fixed << parallelTime << setprecision(6) << " sec" << endl; 
    cout << "Max value: " << max2 << " and its index is: " << maxIndex2 << endl << endl;
    cout << "SPEED UP: " << serialTime / parallelTime << endl;

    return 0;
}