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
    float a = 100.0;
    timeval startTime, endTime;
    double serialTime, parallelTime;

    srand((unsigned int)time(NULL));
    for (int i = 0; i < ARRAY_SIZE; i++)
        array[i] = (float(rand())/float((RAND_MAX)) * a);

    //Serial
    float sumSerial[4] = {0.0, 0.0, 0.0, 0.0};
    float varianceSerial[4] = {0.0, 0.0, 0.0, 0.0};

    float stdDeviationSerial = 0;
    float meanSerial = 0;
    gettimeofday(&startTime, NULL);

    for(int i = 0; i < ARRAY_SIZE; i += 4)
        sumSerial[0] += array[i];
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        sumSerial[1] += array[i + 1];
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        sumSerial[2] += array[i + 2];
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        sumSerial[3] += array[i + 3];
    meanSerial = (sumSerial[0] + sumSerial[1] + sumSerial[2] + sumSerial[3]) / ARRAY_SIZE;

    for(int i = 0; i < ARRAY_SIZE; i += 4)
        varianceSerial[0] += ((array[i] - meanSerial) * (array[i] - meanSerial));
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        varianceSerial[1] += ((array[i + 1] - meanSerial) * (array[i + 1] - meanSerial));
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        varianceSerial[2] += ((array[i + 2] - meanSerial) * (array[i + 2] - meanSerial));
    for(int i = 0; i < ARRAY_SIZE; i += 4)
        varianceSerial[3] += ((array[i + 3] - meanSerial) * (array[i + 3] - meanSerial));
    stdDeviationSerial = sqrt((varianceSerial[0] + varianceSerial[1] + varianceSerial[2] + varianceSerial[3]) / ARRAY_SIZE);

    gettimeofday(&endTime, NULL); 
    serialTime = (endTime.tv_sec - startTime.tv_sec) * 1e6; 
    serialTime = (serialTime + (endTime.tv_usec - startTime.tv_usec)) * 1e-6; 

    //Parallel
    __m128 sumParallel = _mm_set_ps(0.0, 0.0, 0.0, 0.0), varianceParall = _mm_set_ps(0.0, 0.0, 0.0, 0.0);
    __m128 values;
    float totalSumParallel = 0, meanParall = 0, totalVarianceParallel = 0;
    float stdDeviationParallel = 0;
    gettimeofday(&startTime, NULL);
    
    for(int i = 0; i < ARRAY_SIZE; i += 4){
        values = _mm_loadu_ps(&array[i]);
        sumParallel = _mm_add_ps(values, sumParallel);
    }
    sumParallel = _mm_hadd_ps(sumParallel, sumParallel);
    sumParallel = _mm_hadd_ps(sumParallel, sumParallel);
    totalSumParallel = _mm_cvtss_f32(sumParallel);
    meanParall = totalSumParallel / ARRAY_SIZE;
    __m128 meanHelper = _mm_set1_ps(meanParall);

    for(int i = 0; i < ARRAY_SIZE; i += 4){
        values = _mm_loadu_ps(&array[i]);
        values = _mm_sub_ps(values, meanHelper);
        values = _mm_mul_ps(values, values);
        varianceParall = _mm_add_ps(varianceParall, values);
    }
    varianceParall = _mm_hadd_ps(varianceParall, varianceParall);
    varianceParall = _mm_hadd_ps(varianceParall, varianceParall);
    totalVarianceParallel = _mm_cvtss_f32(varianceParall);
    stdDeviationParallel = sqrt( totalVarianceParallel / ARRAY_SIZE);

    gettimeofday(&endTime, NULL); 
    parallelTime = (endTime.tv_sec - startTime.tv_sec) * 1e6; 
    parallelTime = (parallelTime + (endTime.tv_usec - startTime.tv_usec)) * 1e-6;

    cout << "Parall Programming - CA2 - Q2 - 810196512 - 810196469" << endl << endl;
    cout << "Serial Time: " << fixed << serialTime << setprecision(6) << " sec" << endl;
    cout << "Mean: " << meanSerial << endl;
    cout << "Standard Deviation: " << stdDeviationSerial << endl << endl; 
    cout << "Parallel Time: " << fixed << parallelTime << setprecision(6) << " sec" << endl;
    cout << "Mean: " << meanParall << endl;
    cout << "Standard Deviation: " << stdDeviationParallel << endl << endl; 
    cout << "SPEED UP: " << serialTime / parallelTime << endl;

    return 0;
}