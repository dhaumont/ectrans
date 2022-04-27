#define rocfftSafeCall(err) __rocfftSafeCall(err, __FILE__, __LINE__)
#include <rocfft.h>
#include <stdio.h>

static const char *_rocfftGetErrorEnum(rocfft_status error)
{
  switch (error)
    {
#define rs(x) case rocfft_status_##x: return #x
	rs(success);
	rs(failure);
	rs(invalid_arg_value);
	rs(invalid_dimensions);
	rs(invalid_array_type);
	rs(invalid_strides);
	rs(invalid_distance);
	rs(invalid_offset);
	rs(invalid_work_buffer);
#undef rs
    }
  return "unknown";
}


inline void __rocfftSafeCall(rocfft_status err, const char *file, const int line)
{
  if (rocfft_status_success != err) 
    {
      fprintf(stderr, "rocfft at '%s:%d'\n", file, line);
      fprintf(stderr, "rocfft error %d: %s\nterminating!\n", err, _rocfftGetErrorEnum (err));
	  fflush(stderr);
      // cudaDeviceReset(); not sure what this means for AMD device
    }
}

extern "C"
void
#ifdef TRANS_SINGLE
execute_plan_fftr_c_(rocfft_plan *PLANp, float *data_in, float *data_out)
#else
execute_plan_fftr_c_(rocfft_plan *PLANp, double *data_in, double *data_out)
#endif
{
	
// #ifdef TRANS_SINGLE
	// fprintf(stderr, "\n\n(single precision)\n");
// #else
	// fprintf(stderr, "\n\n(double precision)\n");
// #endif
	// fprintf(stderr, "executing plan at %p\n", PLANp);
	// fflush (stderr);

	rocfft_execution_info info = NULL;
	rocfftSafeCall( rocfft_execution_info_create(&info) );
	
	rocfftSafeCall( rocfft_plan_get_print(* PLANp ) );

#pragma omp target data use_device_ptr(data_in,data_out)
{
	rocfftSafeCall( rocfft_execute(*PLANp, // plan
										(void**)&data_in, // in_buffer
										(void**)&data_out, // out_buffer
										info // execution info
									) );
	// fprintf(stderr, "plan at %p executed.\n", PLANp);
	// fflush (stderr);
									
}

}


extern "C"
void
#ifdef TRANS_SINGLE
execute_plan_fftr_inplace_c_(rocfft_plan *PLANp, float *data_in)
#else
execute_plan_fftr_inplace_c_(rocfft_plan *PLANp, double *data_in)
#endif
{
	
	
// #ifdef TRANS_SINGLE
	// fprintf(stderr, "\n\n(single precision)\n");
// #else
	// fprintf(stderr, "\n\n(double precision)\n");
// #endif
	// fprintf(stderr, "using plan at %p\n", PLANp);
	// fflush (stderr);

	rocfft_execution_info info = NULL;
	rocfftSafeCall( rocfft_execution_info_create(&info) );
	
	rocfftSafeCall( rocfft_plan_get_print(* PLANp ) );

//#pragma omp target data use_device_ptr(data_in)
{
	rocfftSafeCall( rocfft_execute(*PLANp, // plan
										(void**)&data_in, // in_buffer
										(void**)&data_in, // out_buffer
										info // execution info
									) );

	// fprintf(stderr, "plan at %p executed.\n", PLANp);
	// fflush (stderr);
}

}

