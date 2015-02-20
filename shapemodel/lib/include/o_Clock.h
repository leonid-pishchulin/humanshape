#ifndef __INCLUDED_TIME
#define __INCLUDED_TIME


#include <ostream>

#ifdef WIN32
	#if !(defined(_WINSOCKAPI_) || defined(_WINSOCK_H))
	/* The check above prevents the winsock2 inclusion if winsock.h already was
	   included, since they can't co-exist without problems */
	#include <winsock2.h>
	#endif
	#include <time.h>
#else
	#include <sys/time.h>
#endif

class o_Clock
{
  struct timeval tv;

public:
  /* gets the current time when it is constructed */
  o_Clock(void);

  /* sets the time (with mod wrap) */
  o_Clock(unsigned long int, double);
  
  /* returns the number of microseconds difference */
  double operator-(const o_Clock&) const;

  /* increases the time (by microseconds) */
  o_Clock operator+(const double) const;
  
  bool operator>(const o_Clock&) const;
  bool operator==(const o_Clock&) const;
  bool operator<(const o_Clock&) const;  

  /* gets the current time */
  void update(void);

  friend std::ostream& operator<<(std::ostream &s, const o_Clock &);
};
#ifdef WIN32
/* windows implementation, coded by someone else */
int gettimeofday( struct timeval* tv, void* timezone );
#endif

#endif
