#include "../include/timer.h"

/**
 * return the # of milliseconds elapsed since the initialization
 */

timer::timer()
{
	gettimeofday(&start, 0);
}

void timer::restart()
{
	*this = timer();
}

size_t timer::get_elapsed()
{
	timeval cur;
	gettimeofday(&cur, 0);
	return diff_timeval(cur, start);
}

size_t timer::diff_timeval(const timeval& end, const timeval& start)
{
	size_t start_clock = start.tv_sec * 1000000LL + start.tv_usec;
	size_t end_clock = end.tv_sec * 1000000LL + end.tv_usec;
	return end_clock - start_clock;
}
