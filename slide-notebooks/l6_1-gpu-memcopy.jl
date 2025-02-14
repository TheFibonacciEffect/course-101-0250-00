#src # This is needed to make this run as normal Julia file
using Markdown #src

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
#nb # _Lecture 6_
md"""
# GPU computing and performance assessment
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
### The goal of this lecture 6 is to tackle:

- GPU architecture and kernel programming
- GPU computing and performance assessment
- Unit testing and reference tests in Julia
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
## GPU architecture and kernel programming
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
We'll get started with a brief overview of the Nvidia GPU architecture and how to program it.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
The Nvidia general purpose GPUs we will use in this course can be programmed using the [CUDA language extension](https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html).

CUDA is accessible in Julia via [CUDA.jl](https://cuda.juliagpu.org/stable/), which exposes most of the native CUDA features to the Julia ecosystem.

Note, however, that CUDA.jl does not use `nvcc`, the Nvidia compiler, but compiles like other Julia code just ahead of time with [LLVM](https://llvm.org). 
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
First, let's distinguish among CPU, GPU, hardware, application and CUDA.

What are _**host**_ and _**device**_?

The _**host**_ is the system CPU. The system memory (DRAM) linked to the CPU is the host memory. The GPU is called a _**device**_ and GPU memory is device memory.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
The GPU hardware is composed of Global (DRAM) memory, L2 cache and many streaming multi-processors (SMs). Each SM contains many compute units (called "CUDA cores" by Nvidia), registers, L1 cache (can be repurposed as shared memory depending on the architecture) and read-only memory. 

> The CUDA programming model provides an abstraction of GPU architecture that acts as a bridge between an application and its possible implementation on GPU hardware. [*[ref]*](https://developer.nvidia.com/blog/cuda-refresher-cuda-programming-model/)
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
In the CUDA programming model, `blocks` of `threads` compose the `grid`. In our implementation, we want to map one thread to each finite-difference cell of the 2D Cartesian domain.

The figure hereafter depicts the relation between the CUDA domain and the finite-difference domain:
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
![cuda_grid](./figures/l6_cuda_grid.png)
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Indices `ix` and `iy` replace the loop indices providing a "vectorised" map of threads - the core to leverage GPU performance. We'll come back to this in a second part of this lecture.

In the CUDA programming model, `blocks` (red) of `threads` compose the `grid` (green).

In our implementations, we will map one thread (red box) to each cell of the 2D Cartesian domain (blue). Other mappings are possible, of course.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
*How does it relate to the GPU hardware?*
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
All threads of a block are guaranteed to be executed concurrently on an SM (yellow box) and therefore share SM resources such as registers, L1 cache (/shared memory) and read-only memory. 

We'll see later that the performance of a GPU application is highly sensitive to the optimal choice of the thread, block, grid layout, the so-called kernel launch parameters.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Writing a Julia GPU function (aka kernel) copying array `A` to array `B` with the layout from the above figure looks as follow
"""
using CUDA

function copy!(A, B)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    A[ix,iy] = B[ix,iy]
    return
end

threads = (4, 3)
blocks  = (2, 2)
nx, ny  = threads[1]*blocks[1], threads[2]*blocks[2]
A       = CUDA.zeros(Float64, nx, ny)
B       =  CUDA.rand(Float64, nx, ny)

@cuda blocks=blocks threads=threads copy!(A, B)
synchronize()

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
_**Playing with GPUs: the rules**_
- Current GPUs allow typically a maximum of 1024 threads per block.
- The maximum number of blocks allowed is huge; computing the largest possible array on the GPU will make you run out of device memory (currently 16-80 GB) before hitting the maximal number of blocks when selecting sensible kernel launch parameters (usually threads per block >= 128).
- Threads, blocks and grid have 3D "Cartesian" topology, which is very useful for 1D, 2D and 3D Cartesian finite-difference domains.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
With this short overview we should have the important concepts in mind to get started with GPU computing 🚀
"""

#nb # > 💡 note: A more complete introduction to CUDA (or refresher) can be accessed [here](https://developer.nvidia.com/blog/tag/cuda-refresher/). Julia GPU resources can be accessed at [https://juliagpu.org](https://juliagpu.org).
#md # \note{A more complete introduction to CUDA (or refresher) can be accessed [here](https://developer.nvidia.com/blog/tag/cuda-refresher/). Julia GPU resources can be accessed at [https://juliagpu.org](https://juliagpu.org).}

#src ######################################################################### 
#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
## GPU computing and performance assessment
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
### The goal of this part is to:
1. Learn about:
  - How to establish the peak memory throughput of your GPU
  - GPU array and kernel programming
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
2. Consolidate:
  - The basics of benchmarking
  - How to compute achieved memory throughput


[*This content is distributed under MIT licence. Authors: S. Omlin (CSCS), L. Räss (ETHZ).*](https://github.com/eth-vaw-glaciology/course-101-0250-00/blob/main/LICENSE.md)
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
In order to get started, we need to connect to a machine which has GPU(s).

Let's take some time to get started.

👉 Getting started:
- Fetch your login infos in the `daint_login.md` file within your personal Polybox folder.
- Head to [Software install](/software_install/#gpu_computing_on_piz_daint) for the directions.
- Finally, fetch the [`l6_1-gpu-memcopy.ipynb`](https://github.com/eth-vaw-glaciology/course-101-0250-00/blob/main/slide-notebooks/notebooks/l6_1-gpu-memcopy.ipynb) notebooks for this lecture and upload them to your `scratch` on Piz Daint.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: Values reported in this notebook are for the Nvidia P100 16GB PCIe GPU.
#md # \note{Values reported in this notebook are for the Nvidia P100 16GB PCIe GPU.}


#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We will use the packages `CUDA` and `BenchmarkTools` to create a little performance laboratory:
"""

import Pkg; Pkg.add("BenchmarkTools");
using CUDA
using BenchmarkTools
#src #using IJulia
#src #using Plots

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
### Scientific applications' performance
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
The performance of most scientific applications nowadays is bound by memory access speed (*memory-bound*) rather than by the speed computations can be done (*compute-bound*).

The reason is that current GPUs (and CPUs) can do many more computations in a given amount of time than they can access numbers from main memory.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
This imbalance can be quantified by dividing the computation peak performance [GFLOP/s] by the memory access peak performance [GB/s] and multiplied by the size of a number in Bytes (for simplicity, theoretical peak performance values as specified by the vendors can be used). For example for the Tesla P100 GPU, it is:

$$ \frac{5300 ~\mathrm{[GFlop/s]}}{732 ~\mathrm{[GB/s]}}~×~8 = 58 $$

(here computed with double precision values taken from [the vendor's product specification sheet](https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/tesla-p100/pdf/nvidia-tesla-p100-PCIe-datasheet.pdf)).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
So we can do 58 floating point operations per number read from main memory or written to it.

As a consequence, we can consider **floating point operations be "for free"** when we work in the memory-bounded regime as in this lecture.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Therefore, let us start with investigating the performance of different ways to express and launch GPU memory copies. We will wrap all of these memory copies in functions, to enable the Julia compiler to optimize them best.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
There exists already the function `copyto!`, which permits to copy data from one pre-allocated array to another; thus, we start with analysing this function's performance.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
But first, let us list what GPUs are available and make sure we assign no more than one user per GPU:
"""
collect(devices())
device!(0) # select a GPU between 0-7


#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
To this purpose, we allocate two arrays and benchmark the function using `BenchmarkTools`:
"""
nx = ny = 32
A = CUDA.zeros(Float64, nx, ny);
B = CUDA.rand(Float64, nx, ny);
@benchmark begin copyto!($A, $B); synchronize() end

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: Previously defined variables are interpolated with `$` into the benchmarked expression.
#md # \note{Previously defined variables are interpolated with `$` into the benchmarked expression.}

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: If not specified otherwise, `CUDA.zeros(nx, ny)` allocates `Float32`.
#md # \warn{If not specified otherwise, `CUDA.zeros(nx, ny)` allocates `Float32`.}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Time samples resulting from benchmarking as just performed follow normally a right skewed distribution.

For such distribution, the median is the most robust of the commonly used estimators of the central tendency; the minimum is in general also a good estimator as hardware cannot by accident run faster than with the ideal and it is as a result commonly used for reporting performance (for more information on estimators see [here](https://juliaci.github.io/BenchmarkTools.jl/stable/manual/#Which-estimator-should-I-use?)).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Using `@belapsed` instead of `@benchmark`, we directly obtain the minimum of the taken time samples:
"""
t_it = @belapsed begin copyto!($A, $B); synchronize() end

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
Now, we know that it does not take "an awful lot of time". Of course, we do not want to stop here, but figure out how good the achieved performance was.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
To this aim, we compute the *total memory throughput*, `T_tot` [GB/s], which is defined as the volume of the copied data [GB] divided by the time spent [s]:
"""
T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: The factor `2` comes from the fact that the data is read and written (`2` operations).
#md # \note{The factor `2` comes from the fact that the data is read and written (`2` operations).}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Compare now `T_tot` with the known peak memory throughput, `T_peak`, which is found e.g. in scientific or vendor publications (for the Nvidia Tesla P100 GPUs, it is 559 GB/s, according to [this source](https://doi.org/10.1109/P3HPC51967.2020.00006).
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: Achievable peak memory throughput is usually significantly lower than the *theoretical peak bandwidth* announced by the vendor (for the [Tesla P100 GPUs](https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/tesla-p100/pdf/nvidia-tesla-p100-PCIe-datasheet.pdf), the latter is 732 GB/s as noted already earlier).
#md # \note{Achievable peak memory throughput is usually significantly lower than the *theoretical peak bandwidth* announced by the vendor (for the [Tesla P100 GPUs](https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/tesla-p100/pdf/nvidia-tesla-p100-PCIe-datasheet.pdf), the latter is 732 GB/s as noted already earlier).}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
#nb # > 💡 note: Here 1 GB is 1e9 Bytes as in the publication, where the peak memory throughput of the Tesla P100 GPU was obtained from.
#md # \note{Here 1 GB is 1e9 Bytes as in the publication, where the peak memory throughput of the Tesla P100 GPU was obtained from.}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
You have surely found `T_tot` to be orders of magnitude below `T_peak`. This is to be expected when copying a small array.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Let us determine how `T_tot` behaves with increasing array sizes:
"""
array_sizes = []
throughputs = []
for pow = 0:11
    nx = ny = 32*2^pow
    if (3*nx*ny*sizeof(Float64) > CUDA.available_memory()) break; end
    A = CUDA.zeros(Float64, nx, ny);
    B = CUDA.rand(Float64, nx, ny);
    t_it = @belapsed begin copyto!($A, $B); synchronize() end
    T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it
    push!(array_sizes, nx)
    push!(throughputs, T_tot)
    println("(nx=ny=$nx) T_tot = $(T_tot)")
#src     #IJulia.clear_output(true)  # Pass wait=True to wait until new output before clearing
#src     #display(plot(array_sizes, throughputs))
    CUDA.unsafe_free!(A)
    CUDA.unsafe_free!(B)
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
You can observe that the best performance is on pair with `T_peak` or a bit lower (measured 522 GB/s with the Tesla P100 GPU) as `copyto!` is a function that needs to work in all possible cases and it is not specifically optimised for a particular hardware.

Furthermore, we note that best performance is obtained for large arrays (in the order of Gigabytes).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We will use the array size for which we obtained the best result for the remainder of the performance experiments:
"""
T_tot_max, index = findmax(throughputs)
nx = ny = array_sizes[index]
A = CUDA.zeros(Float64, nx, ny);
B = CUDA.rand(Float64, nx, ny);

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
### GPU array programming
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
Let us now create our own memory copy function using GPU *Array Programming* (AP).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We can write a memory copy simply as `A .= B`; and wrap it in a function using Julia's concise notation, it looks as follows:
"""
@inbounds memcopy_AP!(A, B) = (A .= B)

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
#nb # > 💡 note: We use `@inbounds` macro to make sure no array bounds checking is performed, which would slow down significantly. Note, furthermore, that outside of these exercises it can be more convenient not to use the `@inbounds` macro, but to deactivate bounds checking instead globally for high performance runs by calling julia as follows : `julia --check-bounds=no ...`
#md # \note{We use `@inbounds` macro to make sure no array bounds checking is performed, which would slow down significantly. Note, furthermore, that outside of these exercises it can be more convenient not to use the `@inbounds` macro, but to deactivate bounds checking instead globally for high performance runs by calling julia as follows : `julia --check-bounds=no ...`}

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: `A = B` would not do a memcopy, but make `A` an alias of `B`, i.e. make `A` point to the same data in memory as `B`.
#md # \note{`A = B` would not do a memcopy, but make `A` an alias of `B`, i.e. make `A` point to the same data in memory as `B`.}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We also benchmark it and compute `T_tot`:
"""
t_it = @belapsed begin memcopy_AP!($A, $B); synchronize() end
T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it

md"""
The performance you observe might be a little lower than with the `copyto!` function (measured 496 GB/s with the Tesla P100 GPU).
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
The few experiments that we have done together so far have shown you already that performing memory copy with maximal possible performance (T_peak) is not a completely trivial task.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
### GPU kernel programming
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
We will now use GPU *Kernel Programming* (KP) to try to get closer to `T_peak`.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
A memory copy kernel can be written e.g. as follows:
"""
@inbounds function memcopy_KP!(A, B)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    A[ix,iy] = B[ix,iy]
    return nothing
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Then, in order to copy the (entire) array `B` to `A`, we need to launch the kernel such that the above indices `ix` and `iy` map exactly to each array cell.

Therefore, we need to have `blocks[1]*threads[1] == nx` and `blocks[2]*threads[2] == ny`.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We will try first with the simplest possible option using only one thread per block:
"""
threads = (1, 1)
blocks  = (nx, ny)
t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_KP!($A, $B); synchronize() end
T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
`T_tot` is certainly orders of magnitude below `T_peak` with this kernel launch parameters.

We need to take into account that single threads cannot run completely independently, but threads are launched in small groups within a block, called *warps*; a warp consists of 32 threads on current GPUs.

Furthermore, warps should access contiguous memory for best performance.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We therefore retry using 32 threads (one warp) per block as follows:
"""
threads = (32, 1)
blocks  = (nx÷threads[1], ny)
t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_KP!($A, $B); synchronize() end
T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
#nb # > 💡 note: For simplicity, the number of threads was set here explicitly to 32; more future proof would be to retrieve the warp size from the corresponding CUDA attribute by doing: `attribute(device(),CUDA.DEVICE_ATTRIBUTE_WARP_SIZE)`.
#md # \note{For simplicity, the number of threads was set here explicitly to 32; more future proof would be to retrieve the warp size from the corresponding CUDA attribute by doing: `attribute(device(),CUDA.DEVICE_ATTRIBUTE_WARP_SIZE)`.}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
`T_tot` is now probably in the order of magnitude of `T_peak`, yet depending on the used GPU it can be still significantly below (measured 302 GB/s with the Tesla P100 GPU).
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
If `T_tot` is significantly below `T_peak`, then we need to set the numbers of threads per block closer to the maximum the GPU allows.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Let us determine how `T_tot` behaves with an increasing number of threads per blocks:
"""
max_threads  = attribute(device(),CUDA.DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK)
thread_count = []
throughputs  = []
for pow = Int(log2(32)):Int(log2(max_threads))
    threads = (2^pow, 1)
    blocks  = (nx÷threads[1], ny)
    t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_KP!($A, $B); synchronize() end
    T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it
    push!(thread_count, prod(threads))
    push!(throughputs, T_tot)
    println("(threads=$threads) T_tot = $(T_tot)")
#src     #IJulia.clear_output(true)
#src     #display(plot(thread_count, throughputs))
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
You should observe now that beyond a certain minimum number of threads per block (64 with the Tesla P100 GPU), `T_tot` is quite close to `T_peak` (which exact thread/block configuration leads to the best `T_tot` depends on the used GPU architecture).
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
Instead of increasing the number of threads only in the x dimension, we can also do so in the y dimension.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We keep though 32 threads in the x dimension in order to let the warps access contiguous memory:
"""
thread_count = []
throughputs  = []
for pow = 0:Int(log2(max_threads/32))
    threads = (32, 2^pow)
    blocks  = (nx÷threads[1], ny÷threads[2])
    t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_KP!($A, $B); synchronize() end
    T_tot = 2*1/1e9*nx*ny*sizeof(Float64)/t_it
    push!(thread_count, prod(threads))
    push!(throughputs, T_tot)
    println("(threads=$threads) T_tot = $(T_tot)")
#src     #IJulia.clear_output(true)
#src     #display(plot(thread_count, throughputs))
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
`T_tot` is even slightly better in general. Much more important is though that a thread block accesses now not a 1D-line of the arrays, but a 2D block.

We will see later that this is of great benefit when, e.g., computing finite difference derivatives in x and y direction.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
So far, we experimented with memory copy in the strict sense: copy an array from one place to the other. When doing computations, we often read more data than we write.

We will therefore also do a few experiments on another commonly benchmarked case: read two arrays and write only one.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
We modify therefore the previous kernel to take a third array `C` as input and add it to `B` (the rest is identical):
"""
@inbounds function memcopy2_KP!(A, B, C)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    A[ix,iy] = B[ix,iy] + C[ix,iy]
    return nothing
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Then, we test exactly as for the previous kernel how `T_tot` behaves with an increasing number of threads per blocks in y dimension, keeping it fixed to 32 in x dimension:
"""
C = CUDA.rand(Float64, nx, ny);
thread_count = []
throughputs  = []
for pow = 0:Int(log2(max_threads/32))
    threads = (32, 2^pow)
    blocks  = (nx÷threads[1], ny÷threads[2])
    t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy2_KP!($A, $B, $C); synchronize() end
    T_tot = 3*1/1e9*nx*ny*sizeof(Float64)/t_it
    push!(thread_count, prod(threads))
    push!(throughputs, T_tot)
    println("(threads=$threads) T_tot = $(T_tot)")
#src     #IJulia.clear_output(true)
#src     #display(plot(thread_count, throughputs))
end

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
#nb # > 💡 note: There is now a factor `3` instead of `2` in the computation of `T_tot`: `2` arrays are read and `1` written (`3` operations).
#md # \note{There is now a factor `3` instead of `2` in the computation of `T_tot`: `2` arrays are read and `1` written (`3` operations).}

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Compare now the best measured `T_tot` to the `T_peak` obtained from the publication and if it is higher, then it means we need to correct `T_peak` to take the value of the `T_tot` measured (`T_tot` measured with the Tesla P100 GPU is 561 GB/s, i.e., 2 GB/s higher than the `T_peak` obtained from the publication mentioned earlier).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Note that the `T_peak` reported in the publication was obtained with a slightly different kernel which multiplies C with a scalar in addition; it is usually referred to as *triad*.

For completeness, we will also quickly benchmark a *triad* kernel.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
To this purpose, we will directly use the best thread/block configuration that we have found in the previous experiment:
"""
@inbounds function memcopy_triad_KP!(A, B, C, s)
    ix = (blockIdx().x-1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y-1) * blockDim().y + threadIdx().y
    A[ix,iy] = B[ix,iy] + s*C[ix,iy]
    return nothing
end

s = rand()

T_tot_max, index = findmax(throughputs)
threads = (32, thread_count[index]÷32)
blocks  = (nx÷threads[1], ny÷threads[2])
t_it = @belapsed begin @cuda blocks=$blocks threads=$threads memcopy_triad_KP!($A, $B, $C, $s); synchronize() end
T_tot = 3*1/1e9*nx*ny*sizeof(Float64)/t_it

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
There should be no significant difference between `T_tot` of this triad kernel and of the previous kernel (with the Tesla P100 GPU, it is 561 GB/s with both kernels).
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Finally, let us also check the triad performance we obtain with GPU array programming:
"""
@inbounds memcopy_triad_AP!(A, B, C, s) = (A .= B.+ s.*C)

t_it = @belapsed begin memcopy_triad_AP!($A, $B, $C, $s); synchronize() end
T_tot = 3*1/1e9*nx*ny*sizeof(Float64)/t_it

#nb # > 💡 Note that in this memory copy function we use everywhere dots to fuse vectorized operations and avoid any allocations of temporary arrays (see [here](https://docs.julialang.org/en/v1/manual/performance-tips/#More-dots:-Fuse-vectorized-operations) for more information).

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
`T_tot` is probably a bit lower than in the above experiment, but still rather close to `T_peak`.
"""

#src ######################################################################### 
#nb # %% A slide [markdown] {"slideshow": {"slide_type": "slide"}}
md"""
Congratulations! You have successfully made it through the memory copy kernel optimization experiments and learned about the fundamental parameters determining memory throughput.
"""

#nb # %% A slide [markdown] {"slideshow": {"slide_type": "fragment"}}
md"""
One moment! For the following exercises you will need the parameters we have established here for best memory access:
"""
println("nx=ny=$nx; threads=$threads; blocks=$blocks")

