inline int gpuGetMaxGflopsDeviceId()
{
  int current_device = 0;
  int sm_per_multiproc = 0;
  int max_perf_device = 0;
  int device_count = 0;
  int best_SM_arch = 0;
  int devices_prohobited = 0;

  unsigned long long max_compute_perf = 0;
  cudaDeviceProp deviceProp;

  cudaGetDeviceCount( &device_count );
  if( device_count == 0 ) {
    printf("gpuGetMaxGflopsDeviceId() ERROR: no devices supporting CUDA.\n");
    exit( EXIT_FAILURE );
  }
}


inline int findCudaDevice( int argc, char **argv )
{
  cudaDeviceProp deviceProp;
  int devID = 0;

  // pick the device with highest Gflops/s
  devID = gpuGetMaxGflopsDeviceId();
}
