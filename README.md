# SpiralReco

Non-Cartesian reconstruction
	Gridding
	NUFFT
	Iterative with regularized
	Exact
Parallel imaging
	coil combination:
		Coil sensitivity map calculation
		Adaptive combine
		SOS
		SOS optimized (noise decorrelation matrix)
	Under sampling
		CG-SENSE
		Coil Compression
		SPIRiT
		ESPIRiT
		g-factor maps
Trajectory correction
	Retrospective correction
		Rotating gradients to physical axis
		GIRF from triangle pulses
		GIRF from  chirp pulses
		Effects of Pre-phasor gradients
	Concurrent correction
		Undo Siemens B0 drift compensation
		Synchronization
	Concomitant field correction
	Density compensation 
		Jackson DCF
		Voronoi.
B0 correction
	Field map calculation
		3-ECHO unwrapping
		UMPIRE
		Regularized
	B_0  correction techniques
		Multi-frequency interpolation
		Semi-Automatic methods
		Iterative (CG-SENSE)
Validation tools
	Reconstruction metric
		SSIM
		PSNR
		TSNR
Distortion maps
