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
destroy_plan_fftr_c_(rocfft_plan *PLANp)
{

rocfftSafeCall( rocfft_plan_destroy(*PLANp) );


}

