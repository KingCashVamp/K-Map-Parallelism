struct Lock {
    int* mutex;

    Lock() {
        int state = 0;
        cudaMalloc((void**)&mutex, sizeof(int));
        cudaMemset(mutex, 0, sizeof(int)); // Initialize mutex to 0
        cudaMemcpy(mutex, &state, sizeof(int), cudaMemcpyHostToDevice);
    }

    ~Lock() {
        cudaFree(mutex);
    }

    __device__ void lock() {
        while (atomicCAS(mutex, 0, 1) != 0) {
            // Busy-wait until lock is acquired
        }
    }

    __device__ void unlock() {
        atomicExch(mutex, 0);
    }
};