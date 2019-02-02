#include <stdio.h>
#include <unistd.h>
#include <linux/unistd.h>
#include <sys/syscall.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sched.h>
#include <signal.h>
#include <setjmp.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include "gt_include.h"


#define ROWS 256
#define COLS ROWS
#define SIZE COLS

#define NUM_CPUS 4
#define NUM_GROUPS NUM_CPUS
#define PER_GROUP_COLS (SIZE/NUM_GROUPS)

#define NUM_THREADS 128
#define PER_THREAD_ROWS (SIZE/NUM_THREADS)

int choice;
static int temp = 1;
long int runtime[NUM_THREADS];
long int exetime[NUM_THREADS];
long int mean_runtime[NUM_THREADS/8];
long int mean_exetime[NUM_THREADS/8];

double stddev_runtime[NUM_THREADS];
double stddev_exetime[NUM_THREADS];

extern void gt_yield();
/* A[SIZE][SIZE] X B[SIZE][SIZE] = C[SIZE][SIZE]
 * Let T(g, t) be thread 't' in group 'g'. 
 * T(g, t) is responsible for multiplication : 
 * A(rows)[(t-1)*SIZE -> (t*SIZE - 1)] X B(cols)[(g-1)*SIZE -> (g*SIZE - 1)] */

typedef struct matrix
{
	int m[SIZE][SIZE];

	int rows;
	int cols;
	unsigned int reserved[2];
} matrix_t;


typedef struct __uthread_arg
{
	matrix_t *_A, *_B, *_C;
	unsigned int reserved0;
        
        int credits; //credit scheduler
        int size;    //credit scheduler

	unsigned int tid;
	unsigned int gid;
	int start_row; /* start_row -> (start_row + PER_THREAD_ROWS) */
	int start_col; /* start_col -> (start_col + PER_GROUP_COLS) */
 	
}uthread_arg_t;
	
struct timeval tv1;

static void generate_matrix(matrix_t *mat, int val,int mat_size)
{

	int i,j;
	mat->rows = mat_size;
	mat->cols = mat_size;
	for(i = 0; i < mat->rows;i++)
		for( j = 0; j < mat->cols; j++ )
		{
			mat->m[i][j] = val;
		}
	return;
}

static void print_matrix(matrix_t *mat, int mat_size)
{
	int i, j;

	for(i=0;i<mat_size;i++)
	{
		for(j=0;j<mat_size;j++)
			printf(" %d ",mat->m[i][j]);
		printf("\n");
	}

	return;
}

static void * uthread_mulmat(void *p)
{
	int i, j, k;
	int start_row, end_row;
	int start_col, end_col;
	unsigned int cpuid;
	struct timeval tv2;
	struct timeval indiv_runtime;
#define ptr ((uthread_arg_t *)p)

	i=0; j= 0; k=0;

	start_row = ptr->start_row;
	end_row = (ptr->start_row + PER_THREAD_ROWS);

#ifdef GT_GROUP_SPLIT
	start_col = ptr->start_col;
	end_col = (ptr->start_col + PER_THREAD_ROWS);
#else
	start_col = 0;
	end_col = SIZE;
#endif

#ifdef GT_THREADS
	cpuid = kthread_cpu_map[kthread_apic_id()]->cpuid;
	fprintf(stderr, "\nThread(id:%d, group:%d, cpu:%d) started",ptr->tid, ptr->gid, cpuid);
#else
	fprintf(stderr, "\nThread(id:%d, group:%d) started",ptr->tid, ptr->gid);
#endif

	for(i = 0/*start_row*/; i < ptr->size/*end_row*/; i++)
		for(j = 0/*start_col*/; j < ptr->size/*end_col*/; j++)
			for(k = 0; k < ptr->size; k++)
				{ptr->_C->m[i][j] += ptr->_A->m[i][k] * ptr->_B->m[k][j]; if(temp == 1 && ptr->tid == 72) {temp = 0;  gt_yield();}
}//Testing yield for a particular scenario with threadid = 72
#ifdef GT_THREADS
	fprintf(stderr, "\nThread(id:%d, group:%d, cpu:%d) finished (TIME : %lu s and %lu us)",
			ptr->tid, ptr->gid, cpuid, (tv2.tv_sec - tv1.tv_sec), (tv2.tv_usec - tv1.tv_usec));
#else
	gettimeofday(&tv2,NULL);

	timersub(&tv2,&tv1, &indiv_runtime);
	runtime[ptr->tid] = indiv_runtime.tv_sec * 1000000L + indiv_runtime.tv_usec;

	fprintf(stderr, "\nThread(id:%d, group:%d) finished (TIME : %lu s and %lu us)",
			ptr->tid, ptr->gid, indiv_runtime.tv_sec, indiv_runtime.tv_usec);
#endif

#undef ptr
	return 0;
}

matrix_t A, B, C;

static void init_matrices(int mat_size)
{
	generate_matrix(&A, 1, mat_size);
	generate_matrix(&B, 1, mat_size);
	generate_matrix(&C, 0, mat_size);

	return;
}


uthread_arg_t uargs[NUM_THREADS];
uthread_t utids[NUM_THREADS];

int main(int argc, char *argv[])
{
        choice = atoi(argv[1]);	
	uthread_arg_t *uarg;
	int inx;
        int ii,jj,kk;
	int size_array[4]={32,64,128,256};
	int credits[4]={25,50,75,100};

        gtthread_app_init();

	if(choice == 0)
	{
	
	init_matrices(SIZE);

	gettimeofday(&tv1,NULL);

	for(inx=0; inx<NUM_THREADS; inx++)
	{
		uarg = &uargs[inx];
		uarg->_A = &A;
		uarg->_B = &B;
		uarg->_C = &C;

		uarg->tid = inx;
		uarg->size = SIZE;
		uarg->gid = (inx % NUM_GROUPS);

		uarg->start_row = (inx * PER_THREAD_ROWS);
#ifdef GT_GROUP_SPLIT
		/* Wanted to split the columns by groups !!! */
		uarg->start_col = (uarg->gid * PER_GROUP_COLS);
#endif

		uthread_create(&utids[inx], uthread_mulmat, uarg, uarg->gid);
	}
	}
	else if(choice == 1)
	{
	
	gettimeofday(&tv1,NULL);


	for(ii=0;ii<4;ii++) //matrix_size
	{
                if(ii == 0)
                {
			init_matrices(32);
                }
		else if(ii == 1)
		{
			init_matrices(64);
		}
		else if(ii == 2)
		{
			init_matrices(128);
		}
		else if(ii == 3)
		{
			init_matrices(256);
		}

		for(jj=0;jj<4;jj++) //credits
		{
		for(kk=0;kk<8;kk++)   //4*4*8 = 128
		{
		inx = 32*ii+8*jj+kk;
		uarg = &uargs[inx];
		uarg->_A = &A;
		uarg->_B = &B;
		uarg->_C = &C;
		
		uarg->credits = credits[jj];
		uarg->size = size_array[ii];

		uarg->tid = inx;

		uarg->gid = (inx % NUM_GROUPS);

		uarg->start_row = (inx * PER_THREAD_ROWS);
#ifdef GT_GROUP_SPLIT
		// Wanted to split the columns by groups !!! 
		uarg->start_col = (uarg->gid * PER_GROUP_COLS);
#endif
		uthread_create(&utids[inx], uthread_mulmat, uarg, uarg->gid, uarg->credits);
	}
	}
	} //3 for loop close


	}
	gtthread_app_exit();

	//print_matrix(&C,256);
	// fprintf(stderr, "********************************");


	if(choice == 1)
	{
	//	for(ii=0;ii<128;ii++)
	//	printf("\nThread:%d  Execution time:%ld  Run time:%ld",ii,exetime[ii],runtime[ii]);	
	//Run time of a thread is from creation to completion
	//Execution time is the cpu time	
	
	
	for(ii=0;ii<4;ii++)
	{
	for(jj=0;jj<4;jj++)
	{
	for(kk=0;kk<8;kk++)
	{
		mean_runtime[ii*4+jj] += runtime[32*ii+8*jj+kk];
		mean_exetime[ii*4+jj] += exetime[32*ii+8*jj+kk];	
	}
	}
	}

	for(ii=0;ii<16;ii++)
	{
		mean_runtime[ii] = mean_runtime[ii]/8;
		mean_exetime[ii] = mean_exetime[ii]/8;
	}
	
	for(ii=0;ii<4;ii++)
	{
	for(jj=0;jj<4;jj++)
	{
	for(kk=0;kk<8;kk++)
	{
		stddev_runtime[ii*4+jj] +=
  pow(abs( (double)runtime[32*ii+8*jj+kk]-(double)mean_runtime[4*ii+jj] ),2);
		stddev_exetime[ii*4+jj] += 
  pow(abs( (double)exetime[32*ii+8*jj+kk]-(double)mean_exetime[4*ii+jj] ),2);

	}
	}
	}

	for(ii=0;ii<16;ii++)
	{
		stddev_runtime[ii] = sqrt(stddev_runtime[ii]/8);
		stddev_exetime[ii] = sqrt(stddev_exetime[ii]/8);
	}
	
	for(ii=0;ii<16;ii++)
	{
	 printf("\nMean Run time for %d is %d Micro seconds",ii,mean_runtime[ii]);
	 printf("\nMean Exe time for %d is %d Micro seconds",ii,mean_exetime[ii]);
	 printf("\nSD Run time for %d is %f Micro seconds",ii,stddev_runtime[ii]);
	 printf("\nSD Exe time for %d is %f Micro seconds",ii,stddev_exetime[ii]);
printf("\n");	
	}

	}//choice end
	return(0);
}
