use messenger

#define FATAL 3
#define NONFATAL 2
#define WARN 1
#define MSG 0

#define signal(ierr,msg) message(ierr,msg,__FILE__,__LINE__)

#define AllocationError  message(FATAL,'Memory allocation error.',__FILE__,__LINE__)
