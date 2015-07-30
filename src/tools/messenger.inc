use messenger

#define FATAL 3
#define NONFATAL 2
#define WARN 1

#define signal(ierr,msg) signal(ierr,msg,__FILE__,__LINE__)

#define AllocationError signal(3,"Allocation Error")