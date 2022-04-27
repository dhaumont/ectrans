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


static int allocatedWorkspace=0;
static void* planWorkspace;
static int planWorkspaceSize=100*1024*1024; //100MB
 
extern "C"
void
create_plan_fftr_c_(rocfft_plan * *PLANp, int *ISIGNp, int *Np, int *LOTp, int *ISTRIDEp, int *IDISTp, int *OSTRIDEp, int *ODISTp, int *INPLACEp)
{

	if(0) {
		fprintf(stderr,"Creating rocFFT plan\n");
		fprintf(stderr,"%s %d \n","N=",*Np);
		fprintf(stderr,"%s %d \n","LOT=",*LOTp);
		fprintf(stderr,"%s %d \n","ISIGN=",*ISIGNp);
		fprintf(stderr,"%s %d \n","ISTRIDEp=",*ISTRIDEp);
		fprintf(stderr,"%s %d \n","IDISTp=",*IDISTp);
		fprintf(stderr,"%s %d \n","OSTRIDEp=",*OSTRIDEp);
		fprintf(stderr,"%s %d \n","ODISTp=",*ODISTp);
		fprintf(stderr,"%s %d \n","INPLACEp=",*INPLACEp);
		fflush (stderr);
	}

	int ISIGN = *ISIGNp;
	size_t N = *Np;
	int LOT = *LOTp;

	// plan description (necessary for strides etcc)
	rocfft_plan_description plan_description;
	rocfftSafeCall( rocfft_plan_description_create(&plan_description) );
	
	size_t in_offsets[1], out_offsets[1], in_strides[1], out_strides[1];
	in_offsets[0]=0; out_offsets[0]=0;
	in_strides[0]=*ISTRIDEp; out_strides[0]=*OSTRIDEp;
	size_t in_distance, out_distance;
	in_distance=*IDISTp;
	out_distance=*ODISTp;
	
	// fprintf(stderr,"Make sure out_distance (%zu) is NX/2+1 !!!!\n",out_distance);
	
	rocfft_array_type intype = ISIGN < 0 ? rocfft_array_type_real : rocfft_array_type_hermitian_interleaved;
	rocfft_array_type outtype = ISIGN < 0 ? rocfft_array_type_hermitian_interleaved : rocfft_array_type_real;

	
	rocfftSafeCall( rocfft_plan_description_set_data_layout(
				plan_description,                        // rocfft_plan_description
				intype,                                  // rocfft_array_type in_array_type
				outtype,                                 // rocfft_array_type out_array_type
				in_offsets,                              // const size_t *in_offsets
				out_offsets,                             // const size_t *out_offsets
				1,                                       // size_t in_strides_size
				in_strides,                              // const size_t *in_strides
				in_distance,                             // size_t in_distance
				1,                                       // size_t out_strides_size
				out_strides,                             // const size_t *out_strides
				out_distance                             // size_t out_distance
			) );

	// plan
	*PLANp = new rocfft_plan;
	rocfft_transform_type transform_type = ISIGN < 0 ? rocfft_transform_type_real_forward : rocfft_transform_type_real_inverse;
	
	#ifdef TRANS_SINGLE
	rocfft_precision prec = rocfft_precision_single;
	#else
	rocfft_precision prec = rocfft_precision_double;
	#endif
	
	rocfft_result_placement plac = ( *INPLACEp == 1 ) ? rocfft_placement_inplace : rocfft_placement_notinplace;
	
	rocfftSafeCall( rocfft_plan_create(*PLANp,
								plac,
								transform_type,
								prec,
								1, // Number of dimensions
								&N, // lengths along different dimensions
								LOT, // Number of transforms
								plan_description)  // Description
				);


	rocfftSafeCall( rocfft_plan_description_destroy(plan_description) );
	
	// fprintf(stderr, "\ncreated plan at %p\n\n", *PLANp);
	// fflush (stderr);

}

