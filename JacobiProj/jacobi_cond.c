/*
 XCode MBP Intel
 All complete
 Compile in workbench by gcc jacobi_cond.c -pthread -lm -o executableName
 */


#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>              /* For timing */
#include <sys/time.h>            /* For timing */
#include <sys/resource.h>
#include <string.h>
#include <pthread.h>
#include <stdbool.h>
#include <semaphore.h>


/****************Global****************************/

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

char *filename;                  /* File name of output file */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 1000000;          /* Maximum iterations, a bound to avoid infinite loop */
double final_diff;               /* Temperature difference between iterations at the end */

/* Thread count */
int thr_count = 2;

/* shared variables between threads */
/*************************************************************/
double** u;                   /* Previous temperatures */
double** w;                   /* New temperatures */

int rowsPerThread, rowsPerLastThread, totalCompletedWorkers = 0, whoPrints = 0;
bool masterCancel = false; //flag for master's cancelation of worker threads
double diffInThr[100]; //store diff calculated by each respective threads

pthread_t ptid[100];
pthread_cond_t waitForMaster = PTHREAD_COND_INITIALIZER; //to wait for master's signal
pthread_cond_t waitForWorkers = PTHREAD_COND_INITIALIZER; //to wait for workers' completion
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; //waitForMaster

/**************************************************************/

int main (int argc, char *argv[])
{
   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   struct rusage usage;

   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state (void);

   /* For convenience of other problem size testing */
   if ((argc == 1) || (argc == 4)) {
      if (argc == 4) {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
         thr_count = atoi(argv[3]);
      } // Otherwise use default grid and thread size
   } else {
     printf("Usage: %s [ <rows> <cols> <threads ]>\n", argv[0]);
     exit(-1);
   }

   printf("Problem size: M=%d, N=%d\nThread count: T=%d\n", M, N, thr_count);

   /* Create the output file */
   filename = argv[0];
   sprintf(filename, "%s.dat", filename);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   gettimeofday (&stime, NULL);
   its = find_steady_state();
   gettimeofday (&etime, NULL);

   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   printf("Converged after %d iterations with error: %8.6f.\n", its, final_diff);
   printf("Elapsed time = %8.4f sec.\n", elapsed);

   getrusage(RUSAGE_SELF, &usage);
   printf("Program completed - user: %.4f s, system: %.4f s\n",
      (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0),
    (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
   printf("no. of context switches: vol %ld, invol %ld\n\n",
            usage.ru_nvcsw, usage.ru_nivcsw);

   print_solution (filename, w);
}

/* Allocate two-dimensional array. */
void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;
   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

/* Set initial and boundary conditions */
void initialize_array (double ***u)
{
   int i, j;

   /* Set initial values and boundary conditions */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         (*u)[i][j] = 25.0;      /* Room temperature */
      (*u)[i][0] = 0.0;
      (*u)[i][N-1] = 0.0;
   }

   for (j = 0; j < N; j++) {
      (*u)[0][j] = 0.0;
      (*u)[M-1][j] = 1000.0;     /* Heat source */
   }
}

/* Print solution to standard output or a file */
void print_solution (char *filename, double **u)
{
   int i, j;
   char sep;
   FILE *outfile;

   if (!filename) { /* if no filename specified, print on screen */
      sep = '\t';   /* tab added for easier view */
      outfile = stdout;
   } else {
      sep = '\n';   /* for gnuplot format */
      outfile = fopen(filename,"w");
      if (outfile == NULL) {
         printf("Can't open output file.");
         exit(-1);
      }
   }

   /* Print the solution array */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         fprintf (outfile, "%6.2f%c", u[i][j], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);

}

void swap (double ***a, double ***b){
   double **temp = *a;
   *a = *b;
   *b = temp;
}

/* Entry function of the worker threads */
void *thr_func(void *arg) {

// Worker's logic here
    struct timeval stime, etime;
    struct rusage usage;
    gettimeofday (&stime, NULL);
    
    int its; //Iteration count
    int i, j;
    int start, end;
    int thrNo = (int)arg; //get the thread number
    
    //------------------calculate the start and end rows for the particular thread------------------//
    if (thrNo == 0){
        start = 1;
    }else{
        start = (thrNo*rowsPerThread) + 1;
    }
    
    if ( ((M-2)%thr_count != 0) && (thrNo == thr_count-1)){
        end = (start + rowsPerLastThread) ; //if M-2 is odd
    } else {
        end = (start + rowsPerThread) ; //if M-2 is even
    }
    //------------------calculate the start and end rows for the particular thread------------------//
    
    
    //---------------------------iteration/main logic of the thread/worker---------------------------//
    for (its = 1; its <= max_its; its++) {
        if (!masterCancel){
           for (i = start; i < end; i++) {
              for (j = 1; j < N-1; j++) {
                  w[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                  if ( fabs((w)[i][j] - (u)[i][j]) > diffInThr[thrNo]){
                     diffInThr[thrNo] = fabs((w)[i][j] - (u)[i][j]);
                  }
              }
           }
            
            //--mutex and cond to update totalCompletedWorkers and wait for master's signal--//
            pthread_mutex_lock(&lock);
            
            totalCompletedWorkers++;
            if (totalCompletedWorkers == thr_count){
                pthread_cond_signal(&waitForWorkers);
            }
            
            pthread_cond_wait(&waitForMaster, &lock);
            
            pthread_mutex_unlock(&lock);
            //--mutex and cond to update totalCompletedWorkers and wait for master's signal--//
            
        }else{
            break;
        }
    }
    //---------------------------iteration/main logic of the thread/worker---------------------------//
    
    
    //------------------------------------if cancel from master------------------------------------//
    gettimeofday (&etime, NULL);
    getrusage(RUSAGE_SELF, &usage);
    
    //---------print worker thread statistics---------//
    while (true){
        if (thrNo == whoPrints){ //so will print in chronological order
            printf("Thread %d has completed - user: %.4f s, system: %.4f s\n", thrNo, (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0), (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
            whoPrints++;
            break;
        }
    }
    //---------print worker thread statistics---------//
    
    //------------------------------------if cancel from master------------------------------------//
    
    return 0;
}

double getThrMaxDiff();
void thrMaxDiffClear();

int find_steady_state (void)
{

// Thread creation and the main control logic
    struct timeval stime, etime;
    struct rusage usage;
    gettimeofday (&stime, NULL);
    
    
    int its = 0;
    double diff = 0.0;
    int iret;
    
    //--------------------------calculate how many rows per thread---------------------------//
    if ( (M-2)%thr_count != 0){ //if rows cannot be divided equally amongst threads
        rowsPerThread = floor((double)(M-2)/thr_count);
        rowsPerLastThread = (M-2) - (rowsPerThread*(thr_count-1));
    } else {
        rowsPerThread = (M-2) /thr_count;
    }
    //---------------------------calculate how many rows per thread---------------------------//
    
    
    //------------------------------------thread creation------------------------------------//
    for (int threadNo=0; threadNo<thr_count; threadNo++){
        ptid[threadNo]= '0' + threadNo;
        iret = pthread_create(&ptid[threadNo], NULL, &thr_func, threadNo);
    }
    //------------------------------------thread creation------------------------------------//
    
    
    //------------------------------------master logic------------------------------------//
    while (its < max_its){
        diff = 0.0;
        
        pthread_mutex_lock(&lock);
        pthread_cond_wait(&waitForWorkers, &lock); //wait for all workers to complete processing
        
        diff = getThrMaxDiff(); //to get the max from diffInThr[]
        thrMaxDiffClear(); //to clear the array diffInThr[]
        
        swap(&u, &w);
        
        if (diff <= EPSILON){
            masterCancel = true;
            pthread_cond_broadcast(&waitForMaster);
            totalCompletedWorkers = 0;
            pthread_mutex_unlock(&lock);
            break;
        } else {
            pthread_cond_broadcast(&waitForMaster);
            totalCompletedWorkers = 0;
            pthread_mutex_unlock(&lock);
        }
        ++its;
    }
    //------------------------------------master logic------------------------------------//
    
    
    //------------------------------------thread joining------------------------------------//
    for (int threadNo=0; threadNo<thr_count; threadNo++){
        pthread_join(ptid[threadNo], NULL);
    }
    //------------------------------------thread joining------------------------------------//
    
    
    gettimeofday (&etime, NULL);
    getrusage(RUSAGE_SELF, &usage);
    
    //---------------------------print master thread statistics---------------------------//
    printf("find_steady_state - user: %.4f s, system: %.4f s\n", (usage.ru_utime.tv_sec + usage.ru_utime.tv_usec/1000000.0), (usage.ru_stime.tv_sec + usage.ru_stime.tv_usec/1000000.0));
    //---------------------------print master thread statistics---------------------------//
    
    final_diff = diff;
    return its;
}

double getThrMaxDiff(){
    double max = 0.0;
    for (int i = 0; i < thr_count; i++)
            if (diffInThr[i] > max)
                max = diffInThr[i];
    return max;
}

void thrMaxDiffClear(){
    for (int i = 0; i < thr_count; i++)
        diffInThr[i] = 0.0;
}
