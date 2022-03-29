#ifndef TIME_STATISTICS_H
#define TIME_STATISTICS_H

// structure used for storing test time measurements
struct time_statistics
{
	float memcpy_h_to_d;
	float memcpy_d_to_h;
	float cuda_algorithm;
	float cpu_algorithm;
	float data_read;
	float data_print;
	int n;
	int l;
};

#endif