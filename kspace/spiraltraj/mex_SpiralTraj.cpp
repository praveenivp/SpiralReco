#include "mex.hpp"
#include "mexAdapter.hpp"
#include "MatlabDataArray.hpp"
#include <iostream>
#include<vector>

#include "vdspiral.h"


class MexFunction : public matlab::mex::Function {


public:

	void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
		//         checkArguments(outputs, inputs);
	
        matlab::data::ArrayFactory factory;
		//get all parameters for vdspiral design

		int Nitlv = (int32_t)inputs[0][0]; 

		double dResolution = inputs[1][0]; //mm

		matlab::data::TypedArray<double> mTA_fov = inputs[2]; //mm
        matlab::data::TypedArray<double> mTA_radius = inputs[3];
        double dGradMaxAmpl = inputs[4][0]; //mT/m
        double dMinRiseTime =inputs[5][0]; //us   //1000. / Slew(mT/m/ms);
		vdspiral::eSpiralType typ = vdspiral::eSpiralType(static_cast<int>(inputs[6][0]));
        double gammabar = 42.5766;
        double grad_raster_time = 10.;
		//peNC_spiral1_fa30_TR100_23052019.edb
        int nfov =4;
		std::vector<double> v_fov({mTA_fov[0],mTA_fov[1],mTA_fov[2],mTA_fov[3]});
        std::vector<double> v_radius({mTA_radius[0],mTA_radius[1],mTA_radius[2],mTA_radius[3]});
		
		


        vdspiral SpiralTraj;
        SpiralTraj.prep(Nitlv,dResolution,v_fov,v_radius,dGradMaxAmpl,dMinRiseTime,typ,gammabar,grad_raster_time);

        std::vector<float> vfGx = SpiralTraj.getGradX(); 
        std::vector<float> vfGy = SpiralTraj.getGradY();    
        double amp_x = SpiralTraj.getMaxAmplitudeX(); //mT/m
        double amp_y = SpiralTraj.getMaxAmplitudeY();//mT/m
        
        matlab::data::TypedArray<double> gradx=factory.createArray<double> ({vfGx.size(),2});
        for (auto it=0;it<vfGx.size();it++)
        {
            gradx[it][0]=amp_x*static_cast<double>(vfGx[it]);//mT/m
            gradx[it][1]=amp_y*static_cast<double>(vfGy[it]);//mT/m
        }

        matlab::data::StructArray st_mom=factory.createStructArray({1},{"PreMomentum","PostMomentum"});
        st_mom[0]["PreMomentum"]=factory.createArray<double>({1,2},{SpiralTraj.getPreMomentumX(),SpiralTraj.getPreMomentumY()});
        st_mom[0]["PostMomentum"]=factory.createArray<double>({1,2},{SpiralTraj.getPostMomentumX(),SpiralTraj.getPostMomentumY()});

        outputs[0]=std::move(gradx);
        outputs[1]=std::move(st_mom);

	

	}
	};
