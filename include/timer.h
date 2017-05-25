#ifndef _H_TIMER
#define _H_TIMER

#include <sys/time.h>
#include <cstdio>

class timer
{
private:
	timeval start;
public:
	timer();
	void restart();
	size_t get_elapsed();
	size_t diff_timeval(const timeval &, const timeval &);
};

#endif /* TIMER_H_ */
