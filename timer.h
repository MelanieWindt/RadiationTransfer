#ifndef __TIMER_H__
#define __TIMER_H__

#include<sys/time.h>

class Timer {
	timeval myTime;
public:
	Timer () {
		gettimeofday (&myTime, 0);
	}

	double stopAndGetElapsedTime () {
		long long startTimeMS = myTime.tv_sec * 1e6 + myTime.tv_usec;
		gettimeofday (&myTime, 0);
		return (myTime.tv_sec * 1e6 + myTime.tv_usec - startTimeMS)/1e6;
	}
};


#endif