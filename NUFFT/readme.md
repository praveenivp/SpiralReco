## StackofSpiral class
`StackofSpiral` class is a flexible interface class for non-cartesian gridding. It can handle for 2D spiral/3D stack-of-spirals data. It can easily switch between GPU and CPU computation. For B0 correction with MFI/MTI, [`StackofSpiralB0`](./../B0Corr/StackofSpiralsB0.m) class which is the extension of this class can be used.


## Setup
It relies on mainly on three matlab/mex functions and its supporting functions. 
1. nufft_init: NUFFT precomputation
2. nufft: forward operator
3. nufft_adj: inverse operator.

CPU functions work out of the box. For GPU, you need to manually copy  mex files from `gpuNUFFT\gpuNUFFT\@gpuNUFFT\private` after compiling the gpuNUFFT package listed below. Alternatevely, you can delete the private folder and add [gpuNUFFT](https://github.com/andyschwarzl/gpuNUFFT ) to path. The binaries in Spiralreco repo are for windows compiled with cuda 10.1.

1. mex_gpuNUFFT_adj_atomic_f
2. mex_gpuNUFFT_forw_atomic_f
3. mex_gpuNUFFT_forw_f
4. mex_gpuNUFFT_precomp_f



### CPU NUFFT
CPU NUFFT files are from [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/) of Fessler.

### gpuNUFFT 
GPU Regridding of arbitrary 3-D/2-D MRI data. Follow https://github.com/andyschwarzl/gpuNUFFT for more details.

- Andreas Schwarzl - andy.schwarzl[at]gmail.com
- Florian Knoll - florian.knoll[at]nyumc.org