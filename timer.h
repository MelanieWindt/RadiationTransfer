#ifndef __TIMER_H__
#define __TIMER_H__

const float CPU_FREQ = 2159; // MHz

inline uint64_t GET_TICK() {
	register unsigned c, d;
	asm volatile("rdtsc" : "=a" (c), "=d" (d));
	return ((uint64_t)c) | (((uint64_t)d) << 32);
}

#define GET_NS() ((uint64_t)(1e3*GET_TICK())/((uint64_t)CPU_FREQ))

class Timer {
	private: 
		static uint64_t prevTime;

	public:
	static void getDelta (const char* str) {
		uint64_t currTime = GET_NS();
		printf("%s \t %llu \n", str, currTime - prevTime );
		prevTime = currTime;
	}

};

uint64_t Timer::prevTime = GET_NS();
#endif