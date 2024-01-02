# SpiralReco
Reconstruction package for model-based reconstruction of accelerated 3D stack-of-spirals data.

### Features
- GIRF correction
- Density compensation function calculation
- easily switch NUFFT algorithms
	- 2D/3D CPU (from [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/))
	- 3D GPU (from [gpuNUFFT](https://github.com/andyschwarzl/gpuNUFFT.git))
- B0 correction with time/frequency segmentation
- Parallel imaging reconstruction
	- CGSENSE
	- SPIRIT (not tested well)
- B0 map preprocessing functions
	- spatial/temporal unwrap
	- iterative denoising of field maps
- Nifti creation

## Installation
For Basic setup
```matlab
cd(userpath)
!git clone --depth=1 https://github.com/pehses/mapVBVD.git
!git clone --depth=1 https://github.com/praveenivp/SpiralReco.git
addpath(fullfile(pwd,'mapVBVD'))
addpath(genpath(fullfile(pwd,'SpiralReco')))
```

The package also can integrate with [gpuNUFFT](https://github.com/andyschwarzl/gpuNUFFT.git). Refer [here](./NUFFT/readme.md) for additional information. Furthermore, [BART](https://github.com/mrirecon/bart.git) was used for coil sensitivity estimation is not included.  

## Usage
Refer to demos
- demo 1 [[matlab livescript](./demo/demo01_phantomdata.mlx)],[[PDF](./demo/demo01_phantomdata.pdf)]
	- Demonstrates features of this package with 3D phantom data. 
- demo 2 [[matlab livescript](./demo/demo02_Simulation2D.mlx)],[[PDF](./demo/demo02_Simulation2D.pdf)]
	- Demonstrates features of this package with simulated data. 

The Documentation I used for HPC reconstruction and post-processing can be found [here](./doc/ClusterProcessing.md).  

## known issues
- GPU NUFFT and undersampling along two direction has problems.
- more to come

## Acknowledgements 
This package has unmodified/modified code from other open-source packages. 
- CPU NUFFT functions are from [MIRT toolbox](https://web.eecs.umich.edu/~fessler/code/) by Jeff Fessler and his group.
- Trajectory calculation are from [spiraltraj](https://github.com/mrphysics-bonn/spiraltraj.git) by Philipp Ehses and Miki Lustig.
- gpuNUFFT functions are from [gpuNUFFT](https://github.com/andyschwarzl/gpuNUFFT.git).
- SPIRIT3D related functions are adapted from  [ESPIRIT package](http://people.eecs.berkeley.edu/~mlustig/Software.html) by Miki Lustig.
- Small parts from [GIRF](https://github.com/MRI-gradient/GIRF) and [ismrm_sunrise_matlab](https://github.com/hansenms/ismrm_sunrise_matlab.git)
- coding styles along with some snippets from Philipp Ehses and Miki Lustig.
- more to come

## References
Please cite the people who developed all these nice methods. Some are listed here.
1. Fessler JA, Sutton BP. Nonuniform fast fourier transforms using min-max interpolation. IEEE Transactions on Signal Processing. 2003;51(2):560-574. doi:10.1109/TSP.2002.807005

2. Sutton BP, Noll DC, Fessler JA. Fast, iterative image reconstruction for MRI in the presence of field inhomogeneities. IEEE Transactions on Medical Imaging. 2003;22(2):178-188. doi:10.1109/TMI.2002.808360

3. Man LC, Pauly JM, Macovski A. Multifrequency interpolation for fast off-resonance correction. Magnetic Resonance in Medicine. 1997;37(5):785-792. doi:10.1002/mrm.1910370523

4. Pruessmann KP, Weiger M, Börnert P, Boesiger P. Advances in sensitivity encoding with arbitrary k-space trajectories. Magnetic Resonance in Medicine. 2001;46(4):638-651. doi:10.1002/mrm.1241

5. Funai AK, Fessler JA, Yeo DTB, Noll DC, Olafsson VT. Regularized field map estimation in MRI. IEEE Transactions on Medical Imaging. 2008;27(10):1484-1494. doi:10.1109/TMI.2008.923956

6. Lustig M, Kim SJ, Pauly JM. A fast method for designing time-optimal gradient waveforms for arbitrary k-space trajectories. IEEE Transactions on Medical Imaging. 2008;27(6):866-873. doi:10.1109/TMI.2008.922699

7. Vannesjo SJ, Graedel NN, Kasper L, et al. Image reconstruction using a gradient impulse response model for trajectory prediction. Magnetic Resonance in Medicine. 2016;76(1):45-58. doi:10.1002/mrm.25841

9. 	Robson PM, Grant AK, Madhuranthakam AJ, Lattanzi R, Sodickson DK, McKenzie CA. Comprehensive quantification of signal-to-noise ratio and g-factor for image-based and k-space-based parallel imaging reconstructions. Magnetic Resonance in Medicine. 











	
