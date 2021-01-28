/*
return the spiral trajectories.
*/

#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include "spline.h"
#include "mtg_functions_pe.h"
#include<vector>
enum eSpiralType {
	SpiralOut = 1,
	SpiralIn = 2,
	DoubleSpiral = 3,
	SpiralInAndOut = 4
};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


//most of this class is from vdspiral.cpp and vdspiral.h of the peNC_spiral sequence. stripped all unnecessary fucntions and added a few get and set fucntions.
class myclass {



protected:
	std::vector<float> m_vfGx, m_vfGy, m_vfGz;
	double m_dAx, m_dAy, m_dAz, m_dAmp, m_dMomX, m_dMomY, m_dMomZ, m_dPreMomX, m_dPreMomY, m_dPreMomZ, m_dPostMomX, m_dPostMomY, m_dPostMomZ;
	double m_dGradRasterTime;
	double m_dLarmorConst;
	eSpiralType m_eSpiralType;

	double m_dResolution, m_dMaxAmplitude, m_dMinRiseTime;
	int m_Nitlv;
	std::vector<double> m_fov, m_radius;



public:
	myclass(void) // Constructor
		: m_dAx(0.)
		, m_dAy(0.)
		, m_dAz(0.)
		, m_dAmp(0.)
		, m_dMomX(0.)
		, m_dMomY(0.)
		, m_dMomZ(0.)
		, m_dPreMomX(0.)
		, m_dPreMomY(0.)
		, m_dPreMomZ(0.)
		, m_dPostMomX(0.)
		, m_dPostMomY(0.)
		, m_dPostMomZ(0.)
		, m_dLarmorConst(0.)
		, m_dGradRasterTime(0)
		, m_eSpiralType(SpiralOut)
		, m_Nitlv(1)
		, m_dResolution(5.) //mm
	{
		//this->m_eSpiralType = SpiralOut;
		//this->m_Nitlv = 8;
		//this->m_dResolution = 2.;//mm
		//this->m_dLarmorConst = 42.575575;//42.5756;
	}


	~myclass(void) // Destructor
	{}

public:
	void setparameters(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T)
	{
		m_Nitlv = Nitlv;
		m_dResolution = res;
		m_dMaxAmplitude = Gmax;
		m_dMinRiseTime = 1000. / Smax;
		m_eSpiralType = spiralType;
		m_dGradRasterTime = T;
		m_fov.clear();
		m_radius.clear();
		for (int i = 0; i < nfov; i++)
		{
			m_fov.push_back(fov_[i]);
			m_radius.push_back(radius_[i]);
		}
		m_dLarmorConst = 42.5756;// 39.9718263;//42.575575;

	}

	void PrintOtherparameters()
	{
		std::cout << this->m_dAx << std::endl;
		std::cout << this->m_dMomX << std::endl;
		std::cout << this->m_dMomY << std::endl;
		std::cout << this->m_dMomZ << std::endl;
		std::cout << this->m_dPreMomX << std::endl;
		std::cout << this->m_dPreMomY << std::endl;
		std::cout << this->m_dPreMomZ << std::endl;
		std::cout << this->m_dPostMomX << std::endl;
		std::cout << this->m_dPostMomY << std::endl;
		std::cout << this->m_dPostMomZ << std::endl;


	}

	void getall(std::vector<double>& vec)
	{
		
		vec.push_back(this->m_dAx);
		vec.push_back(this->m_dAy);
		vec.push_back(this->m_dAz);
					  
		vec.push_back(this->m_dMomX);
		vec.push_back(this->m_dMomY);
		vec.push_back(this->m_dMomZ);
					  
		vec.push_back(this->m_dPreMomX);
		vec.push_back(this->m_dPreMomY);
		vec.push_back(this->m_dPreMomZ);
					 
		vec.push_back(this->m_dPostMomX);
		vec.push_back(this->m_dPostMomY);
		vec.push_back(this->m_dPostMomZ);

		for (int i = 0; i < this->m_vfGx.size(); i++)
		{
			vec.push_back(this->m_vfGx[i]);
		}
		for (int i = 0; i < this->m_vfGy.size(); i++)
		{
			vec.push_back(this->m_vfGy[i]);
		}
		for (int i = 0; i < this->m_vfGz.size(); i++)
		{
			vec.push_back(this->m_vfGz[i]);
		}
	}

	void printGX()
	{

		for (int i = 0; i < this->m_vfGx.size(); i++)
			std::cout << m_vfGx[i] << "  ";
	}

	void printGY()
	{
		for (int i = 0; i < this->m_vfGy.size(); i++)
			std::cout << m_vfGy[i] << "  ";
	}

	bool vdSpiralDesign(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T) {

		long k; // loop index

		if (nfov != nradius) {
			std::cout << "Error: array fov needs to have same size as array radius" << std::endl;
			return false;
		}
		// we need to make copies of fov and radius since we do a unit transformation
		double *fov = new double[nfov];
		double *radius = new double[nradius];

		double kmax = 5. / res;  // kmax = 1/(2*res) BUT: kmax in 1/cm, res in mm

		double fovmax = 0.;
		for (k = 0; k < nfov; ++k) {
			fov[k] = fov_[k] / 10.; // mm->cm
			if (fov[k] > fovmax)
				fovmax = fov[k];
			radius[k] = kmax * radius_[k];
		}


		//double dr = 1./500. * 1./(fovmax/Nitlv);
		double dr = 1. / 100. * 1. / (fovmax / Nitlv); // a little faster
		long   nr = long(kmax / dr) + 1;

		std::vector<double> x, y, z;
		x.resize(nr, 0.);
		y.resize(nr, 0.);
		z.resize(nr, 0.);

		double theta = 0.;
		for (k = 0; k < nr; k++) {
			double r = k * dr;
			double cFoV = fov[nfov - 1];
			for (int l = 0; l < nfov; ++l) {
				if (r < radius[l]) {
					if (l == 0 || l == nfov - 1)
						cFoV = fov[l];
					else {// linearer übergang // linear transition
						double step = (r - radius[l - 1]) / (radius[l] - radius[l - 1]);
						cFoV = step * fov[l] + (1. - step)*fov[l - 1];
					}
					break;
				}
			}
			x[k] = r * cos(theta);
			y[k] = r * sin(theta);
			if (spiralType == DoubleSpiral)
				theta += M_PI * dr*cFoV / Nitlv;
			else
				theta += 2.*M_PI*dr*cFoV / Nitlv;
		}

		delete[] fov; delete[] radius;

		Gmax /= 10.;   // mT/m    -> G/cm
		Smax /= 10.;   // mT/m/ms -> G/cm/ms
		T /= 1000.; // us      -> ms

		int n;
		double g0 = 0.; // to simplify sequence development, our gradient will start at 0.
		double gfin = 0.; // and end at 0.
		double *gx; double *gy; double *gz;
		//clock_t start = clock();
		minTimeGradientRIV(&x[0], &y[0], &z[0], nr, g0, gfin, Gmax, Smax, T, gx, gy, gz, n, -1.5, m_dLarmorConst / 10.);
		//clock_t end = clock();
		//cout << "start = " << start << "  end = " << end << "  end-start = " << end-start << "  CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;

		// determine max gradient amplitudes
		m_dAx = 0.;
		m_dAy = 0.;
		m_dAmp = 0.;
		for (k = 0; k < n; ++k) {
			if (fabs(gx[k]) > m_dAx)
				m_dAx = fabs(gx[k]);
			if (fabs(gy[k]) > m_dAy)
				m_dAy = fabs(gy[k]);
			double dAmp = sqrt(gx[k] * gx[k] + gy[k] * gy[k]);
			if (dAmp > m_dAmp)
				m_dAmp = dAmp;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		///// pehses: start gradient with zero and ramp down the gradient to zero for good at the end /////
		///////////////////////////////////////////////////////////////////////////////////////////////////
		int addRampUp = 1;
		int addRampDn = 1 + (int)(sqrt(gx[n - 1] * gx[n - 1] + gy[n - 1] * gy[n - 1] + gz[n - 1] * gz[n - 1]) / Smax / T);
		if (spiralType == DoubleSpiral) {
			addRampUp = addRampDn;
			m_vfGx.resize(2 * n + addRampUp + addRampDn + 1, 0.);
			m_vfGy.resize(2 * n + addRampUp + addRampDn + 1, 0.);
			// first value of gradients is zero
			// copy & scale gradient to interval -1...+1
			for (k = 0; k < n; ++k) {
				m_vfGx[k + addRampUp] = (float)(gx[n - k - 1] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[k + addRampUp] = (float)(gy[n - k - 1] / (m_dAy > 0 ? m_dAy : 1.));
				m_vfGx[n + k + addRampUp + 1] = (float)(gx[k] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[n + k + addRampUp + 1] = (float)(gy[k] / (m_dAy > 0 ? m_dAy : 1.));
			}
			n = 2 * n + 1;
		}
		else {
			m_vfGx.resize(n + addRampUp + addRampDn, 0.);
			m_vfGy.resize(n + addRampUp + addRampDn, 0.);
			// first value of gradients is zero
			// copy & scale gradient to interval -1...+1
			for (k = 0; k < n; ++k) {
				m_vfGx[k + addRampUp] = (float)(gx[k] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[k + addRampUp] = (float)(gy[k] / (m_dAy > 0 ? m_dAy : 1.));
			}
		}
		delete[] gx; delete[] gy; delete[] gz;

		// linear ramp up from zero (we should already be close)
		for (k = 0; k < addRampUp; ++k) {
			m_vfGx[k] = (float)(m_vfGx[addRampUp] * (double(k) / addRampUp));
			m_vfGy[k] = (float)(m_vfGy[addRampUp] * (double(k) / addRampUp));
		}
		// linear ramp down to zero at the end (we should already be close)
		for (k = 0; k < addRampDn; ++k) {
			m_vfGx[n + addRampUp + k] = (float)(m_vfGx[n + addRampUp - 1] * (1. - (k + 1.) / addRampDn));
			m_vfGy[n + addRampUp + k] = (float)(m_vfGy[n + addRampUp - 1] * (1. - (k + 1.) / addRampDn));
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////

		// G/cm -> mT/m   
		m_dAx *= 10.;
		m_dAy *= 10.;
		m_dAmp *= 10.;

		// now calculate the gradient moments
		m_dMomX = 0.; m_dMomY = 0.; m_dMomZ = 0.;
		for (k = 0; k < (int)m_vfGx.size(); ++k) {
			m_dMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
			m_dMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
		}
		m_dPreMomX = 0.; m_dPreMomY = 0.; m_dPreMomZ = 0.;
		m_dPostMomX = 0.; m_dPostMomY = 0.; m_dPostMomZ = 0.;
		if (spiralType == SpiralIn) {
			m_dPreMomX = m_dMomX;
			m_dPreMomY = m_dMomY;
			m_dPreMomZ = m_dMomZ;
			// we have to time-reverse the trajectory!
			for (long k = 0; k < (int)m_vfGx.size() / 2; ++k) {
				std::swap(m_vfGx[k], m_vfGx[m_vfGx.size() - 1 - k]);
				std::swap(m_vfGy[k], m_vfGy[m_vfGy.size() - 1 - k]);
			}
		}
		else if (spiralType == SpiralOut) {
			m_dPostMomX = m_dMomX;
			m_dPostMomY = m_dMomY;
			m_dPostMomZ = m_dMomZ;
		}
		else if (spiralType == DoubleSpiral) {
			for (k = 0; k < (int)(m_vfGx.size() / 2); ++k) {
				m_dPreMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
				m_dPreMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
			}
			m_dPostMomX = m_dMomX - m_dPreMomX;
			m_dPostMomY = m_dMomY - m_dPreMomY;
			m_dPostMomZ = m_dMomZ - m_dPreMomZ;

		}
		else
			return false;

		return true;
	}


	bool setSpiralType(eSpiralType spiralType) {
		bool bStatus = true;
		if (spiralType != m_eSpiralType) {
			if ((m_eSpiralType != DoubleSpiral) && (spiralType != DoubleSpiral)) {
				// we have to time-reverse the trajectory!
				for (long k = 0; k < (int)m_vfGx.size() / 2; ++k) {
					std::swap(m_vfGx[k], m_vfGx[m_vfGx.size() - 1 - k]);
					std::swap(m_vfGy[k], m_vfGy[m_vfGy.size() - 1 - k]);
				}
				std::swap(m_dPreMomX, m_dPostMomX);
				std::swap(m_dPreMomY, m_dPostMomY);
				std::swap(m_dPreMomZ, m_dPostMomZ);
#ifdef BUILD_SEQU
				bStatus = vdspiral::prepGradients();
#endif
				return bStatus;
			}
			else { // we need to recalculate the trajectory
				m_eSpiralType = spiralType;
				return false; //this->vdSpiralDesign(m_Nitlv, m_dResolution, m_fov, m_radius, m_dMaxAmplitude, m_dMinRiseTime, m_eSpiralType, m_dLarmorConst, m_dGradRasterTime);
			}
		}
		return true;
	}
};




class MexFunction : public matlab::mex::Function {

//protected:
//	//std::vector<float> m_vfGx, m_vfGy, m_vfGz;
//	double m_dAx, m_dAy, m_dAz, m_dAmp, m_dMomX, m_dMomY, m_dMomZ, m_dPreMomX, m_dPreMomY, m_dPreMomZ, m_dPostMomX, m_dPostMomY, m_dPostMomZ;
//	double m_dGradRasterTime;
//	double m_dLarmorConst;

public:

	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		//         checkArguments(outputs, inputs);
		//double multiplier = inputs[0][0];

		//get all parameters for vdspiral design

		int Nitlv = (int32_t)inputs[0][0]; //pMrProt->kSpace().radialViews

		double dResolution = inputs[1][0];

		int nfov = (int32_t)inputs[2][0];
		//peNC_spiral1_fa30_TR100_23052019.edb
		double* fov = new double[nfov];
		double* radius = new double[nfov];
		fov[0] = inputs[3][0]; fov[1] = inputs[3][1]; fov[2] = inputs[3][2]; fov[3] = inputs[3][3];

		radius[0] = inputs[4][0]; radius[1] = inputs[4][1]; radius[2] = inputs[4][2]; radius[3] = inputs[4][3];


		

		double dGradMaxAmpl = inputs[5][0];
		
		
		double Smax = inputs[6][0];  // max slew rate
		double dMinRiseTime = 1000. / Smax;
		eSpiralType typ = eSpiralType((int)inputs[7][0]);
		double GRAD_RASTER_TIME = inputs[8][0];

		std::vector<double> vec_output;


		GetPESpiralTraj(vec_output,Nitlv, dResolution, fov, nfov, radius, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
		//vec_outputs has moments of prewinder, rewinder and actual readout gradients. check getall() for the structure of the vector.


		matlab::data::ArrayFactory a;
		a.createEmptyArray();
		matlab::data::ArrayDimensions dim;
		dim = {2};
		
		matlab::data::TypedArray<double> b = a.createArray<double>({vec_output.size()});

		
		for (int i = 0; i < vec_output.size(); i++)
		{
			b[i] = vec_output[i];
		}


		outputs[0] = std::move(b);
	

	}

	void GetPESpiralTraj(std::vector < double>& out, int Nitlv, double dResolution, double* fov_, int nfov, double* radius_, double& dGradMaxAmpl, double Smax, eSpiralType typ, double GRAD_RASTER_TIME) {

		myclass spiralTraj;
		myclass spiralTraj2;

	/*	out.push_back(fov_[0]);
		out.push_back(fov_[1]);
		out.push_back(fov_[2]);
		out.push_back(fov_[3]);
		out.push_back(radius_[0]);
		out.push_back(radius_[1]);
		out.push_back(radius_[2]);
		out.push_back(radius_[3]);
		out.push_back(Smax);
		out.push_back(GRAD_RASTER_TIME);
		out.push_back((double)typ);*/


		if (typ == SpiralInAndOut)
		{
			// SpiralInAndOut is not convertible to vdspiral::eSpiralType!
			spiralTraj.setparameters(Nitlv, dResolution, fov_, nfov, radius_, nfov, dGradMaxAmpl, Smax, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
			spiralTraj.vdSpiralDesign(Nitlv, dResolution, fov_, nfov, radius_, nfov, dGradMaxAmpl, Smax, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
		}
		else {
			spiralTraj.setparameters(Nitlv, dResolution, fov_, nfov, radius_, nfov, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
			spiralTraj.vdSpiralDesign(Nitlv, dResolution, fov_, nfov, radius_, nfov, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
		}



		if (typ == SpiralInAndOut) {
			// SpiralInAndOut is not convertible to vdspiral::eSpiralType!
			// spiral out
			spiralTraj.setSpiralType(eSpiralType::SpiralIn);

			// spiral in
			spiralTraj2 = spiralTraj;
			spiralTraj.setSpiralType(eSpiralType::SpiralOut);
		}
		else {
			spiralTraj.setSpiralType(typ);
		}


		//out.clear();
		spiralTraj.getall(out);
		if (typ == SpiralInAndOut) 
			spiralTraj2.getall(out);
	}

	/*
	double callvddesign(double& num)
	{

		int n = 2;

		double* fov = new double[n];
		double* rad = new double[n];
		fov[0] = 100;
		fov[1] = 100;
		rad[0] = 50;
		rad[1] = 50;

		eSpiralType typ = SpiralOut;

		double gmax = 30.;
		double T = 1000.;


		vdSpiralDesign(1, 1., fov, n, rad, n, gmax, gmax, typ, T);

		num = num + 1;

		return num;
	}

	*/
	/*
	int main()
	{
		std::cout << "Hello World!\n";


		int Nitlv = 8; //pMrProt->kSpace().radialViews

		double dResolution = 2.;

		int nfov = 4;
		//peNC_spiral1_fa30_TR100_23052019.edb
		double* fov = new double[nfov];
		double* radius = new double[nfov];
		fov[0] = 250; fov[1] = 250; fov[2] = 250; fov[3] = 250;

		radius[0] = 0; radius[1] = 0.15; radius[2] = 0.2; radius[3] = 1;


		eSpiralType typ = eSpiralType::SpiralInAndOut;

		double dGradMaxAmpl = 42.;
		double GRAD_RASTER_TIME = 10.;
		double dMinRiseTime = 5.5;
		double Smax = 1000. / dMinRiseTime;// max slew rate




		myclass spiralTraj;
		myclass spiralTraj2;




		if (typ == SpiralInAndOut)
		{
			// SpiralInAndOut is not convertible to vdspiral::eSpiralType!
			spiralTraj.setparameters(Nitlv, dResolution, fov, nfov, radius, nfov, dGradMaxAmpl, Smax, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
			spiralTraj.vdSpiralDesign(Nitlv, dResolution, fov, nfov, radius, nfov, dGradMaxAmpl, Smax, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
		}
		else {
			//spiralTraj.setparameters(Nitlv, dResolution, fov, nfov, radius, nfov, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
			spiralTraj.vdSpiralDesign(Nitlv, dResolution, fov, nfov, radius, nfov, dGradMaxAmpl, Smax, typ, GRAD_RASTER_TIME);
		}



		if (typ == SpiralInAndOut) {
			// SpiralInAndOut is not convertible to vdspiral::eSpiralType!
			// spiral out
			spiralTraj.setSpiralType(eSpiralType::SpiralIn);

			// spiral in
			spiralTraj2 = spiralTraj;
			spiralTraj.setSpiralType(eSpiralType::SpiralOut);
		}
		else {
			spiralTraj.setSpiralType(typ);
		}



		if (typ != SpiralInAndOut) {
			spiralTraj.printGX();
			std::cout << std::endl << std::endl << std::endl << std::endl;
			spiralTraj.printGY();
		}
		else
		{
			//std::cout << "first part of Spiralinandout traj:" << std::endl;
			std::cout << std::endl << std::endl << std::endl << std::endl;
			spiralTraj.printGX();
			std::cout << std::endl << std::endl << std::endl << std::endl;
			spiralTraj.printGY();
			std::cout << std::endl << std::endl << std::endl << std::endl;
			//std::cout << "Second part of Spiralinandout traj:" << std::endl;
			spiralTraj2.printGX();
			std::cout << std::endl << std::endl << std::endl << std::endl;
			spiralTraj2.printGY();
			std::cout << std::endl << std::endl << std::endl << std::endl;
		}
		spiralTraj.PrintOtherparameters();

	}

	*/

	





	//     void checkArguments(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
	//         std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
	//         matlab::data::ArrayFactory factory;
	// 
	//         if (inputs.size() != 2) {
	//             matlabPtr->feval(u"error", 
	//                 0, std::vector<matlab::data::Array>({ factory.createScalar("Two inputs required") }));
	//         }
	// 
	//         if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE ||
	//             inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE ||
	//             inputs[0].getNumberOfElements() != 1) {
	//             matlabPtr->feval(u"error", 
	//                 0, std::vector<matlab::data::Array>({ factory.createScalar("Input multiplier must be a scalar") }));
	//         }
	// 
	//         if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE ||
	//             inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
	//             matlabPtr->feval(u"error", 
	//                 0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
	//         }
	// 
	//         if (inputs[1].getDimensions().size() != 2) {
	//             matlabPtr->feval(u"error", 
	//                 0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
	//         }
	//     }
};




