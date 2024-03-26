#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>

/*
Paralelizovan i for i while, kao eksponencijalno
*/

double cpu_time(void)
{
  double value;

  value = (double)clock() / (double)CLOCKS_PER_SEC;

  return value;
}

/*
int prime_number(int n)
{
  int i;
  int j;
  int prime;
  int total;

  total = 0;
  short myId, nthreads;
  int start, end, chunk;
#pragma omp parallel default(none)\
                      shared(n, nthreads, chunk)\
                      private(start, end, myId, prime, i, j)\
                      reduction(+:total)
{
  nthreads = omp_get_num_threads();
  myId = omp_get_thread_num();
  if(n<1000)
  { 
    if(myId == 0)
    {
      start = 2;
      end = 2 + n-nthreads;
    }else
    {
      start = end = 2 + n-nthreads + myId;
      if(end>n && start!=n)end=n;
      if(start<2)start = end+1; 
    }
  }else
  {
    chunk = (n-2+nthreads-1) / nthreads;
    start = 2 + myId * chunk;
    end = start+chunk-1 < n ? start+chunk-1 : n;
    if(myId == nthreads-1)end = n;
  }
  for (i = start; i <= end; i++)
  {
    prime = 1;
    for (j = 2; j < i; j++)
    {
      if ((i % j) == 0)
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }
}

  return total;
}
*/

int prime_number(int n)
{
  int i;
  int j;
  int prime;
  int total;

  total = 0;
  #pragma omp parallel default(none) private(i,j, prime) shared(n) reduction(+:total)
  {
    int my_id = omp_get_thread_num();
    int n_threads = omp_get_num_threads();
    
    for (i = 2 + my_id; i <= n; i += n_threads)
    {
      prime = 1;
      for (j = 2; j < i; j++)
      {
        if ((i % j) == 0)
        {
          prime = 0;
          break;
        }
      }
      total = total + prime;
    }
  }
  return total;
}


void timestamp(void)
{
#define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  printf("%s\n", time_buffer);

  return;
#undef TIME_SIZE
}

void test(int n_lo, int n_hi, int n_factor);

int main(int argc, char *argv[])
{
  int n_factor;
  int n_hi;
  int n_lo;

  timestamp();
  printf("\n");
  printf("PRIME TEST\n");

  if (argc != 4)
  {
    n_lo = 1;
    n_hi = 131072;
    n_factor = 2;
  }
  else
  {
    n_lo = atoi(argv[1]);
    n_hi = atoi(argv[2]);
    n_factor = atoi(argv[3]);
  }

  test(n_lo, n_hi, n_factor);

  printf("\n");
  printf("PRIME_TEST\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();

  return 0;
}


// void test(int n_lo, int n_hi, int n_factor)
// {
//   int i;
//   int n;
//   int primes;
//   double ctime;

//   printf("\n");
//   printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
//   printf("\n");
//   printf("         N        Pi          Time\n");
//   printf("\n");

//   n = n_lo;
//   int numRepetitions = 0;
//   int nn=n, nh=n_hi;
//   while(nn<=nh)
//   {
//     numRepetitions++;
//     nn*=n_factor;
//   }

//   int start, end, chunk, myId, nthreads;
// #pragma parallel default(none)\
//                   private(start, end, myId, ctime, primes)\
//                   shared(n, n_hi, n_factor, numRepetitions)
// {
//   myId = omp_get_thread_num();
//   nthreads = omp_get_num_threads();
//   int myReps = 1;

//   if(myId == 0)
//   {
//     myReps = numRepetitions - nthreads + 1;
//   }
//   start = n;
//   for(int i=0; i<myId;i++){
//     start*=n_factor;
//   }
//   for(int j=0; j<myReps; j++)
//   {
//     ctime = cpu_time();

//     primes = prime_number(start);

//     ctime = cpu_time() - ctime;

//     printf("  %8d  %8d  %14f\n", start, primes, ctime);
//     start = start * n_factor;
//   }
// }

//   return; 
// }


void test(int n_lo, int n_hi, int n_factor)
{
  int i;
  int n;
  int primes;
  double ctime;

  printf("\n");
  printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
  printf("\n");
  printf("         N        Pi          Time\n");
  printf("\n");

  n = n_lo;

  while (n <= n_hi)
  {
    ctime = cpu_time();

    primes = prime_number(n);

    ctime = cpu_time() - ctime;

    printf("  %8d  %8d  %14f\n", n, primes, ctime);
    n = n * n_factor;
  }

  return;
}
