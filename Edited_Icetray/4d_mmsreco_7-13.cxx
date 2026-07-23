#include <icetray/I3ConditionalModule.h>
#include <icetray/I3ServiceBase.h>
#include <icetray/I3SingleServiceFactory.h>

#include <dataclasses/I3DOMFunctions.h>
#include <dataclasses/I3Double.h>
#include <icetray/I3Units.h>
#include <dataclasses/geometry/I3Geometry.h>
#include <dataclasses/calibration/I3Calibration.h>
#include <dataclasses/calibration/I3DOMCalibration.h>
#include <dataclasses/physics/I3RecoPulse.h>
#include <simclasses/I3CompressedPhoton.h>
#include <dataclasses/physics/I3Waveform.h>
#include <dataclasses/status/I3DetectorStatus.h>
#include <dataclasses/status/I3DOMStatus.h>
#include <dataclasses/physics/I3Particle.h>


#include <phys-services/I3Calculator.h>
#include <phys-services/I3RandomService.h>

#include <gulliver/I3EventLogLikelihoodBase.h>
#include <gulliver/I3EventHypothesis.h>
#include <cmath>

#include <photospline/splinetable.h>
#include <fstream>
#include <ios>

double C = 299792458; // I3constant

class MMSLikelihood : public I3EventLogLikelihoodBase, public
    I3ServiceBase
{
	public:
		MMSLikelihood(const I3Context &ctx);

                void Configure();
		void SetEvent(const I3Frame &);
		void SetGeometry(const I3Geometry &geo) {};
		double GetLogLikelihood(const I3EventHypothesis &);
		unsigned int GetMultiplicity();

		virtual const std::string GetName() const {
			return I3ServiceBase::GetName();
		}
	private:
		photospline::splinetable<> spline_table_;
		std::string photons_name_;
		bool noise_;
		I3PhotonSeriesMapConstPtr photons_; // Not sure if this is right object
		I3GeometryConstPtr geo_;
};

typedef I3SingleServiceFactory<MMSLikelihood, I3EventLogLikelihoodBase>
    MMSLikelihoodFactory;
I3_SERVICE_FACTORY(MMSLikelihoodFactory);

MMSLikelihood::MMSLikelihood(const I3Context &ctx) : I3ServiceBase(ctx)
{
	AddParameter("SplineTablePath", "Path to spline table", "");
	AddParameter("InputPhotons", "Photons to use", "I3Photons"); // need to ask what arguments are and how specific
	AddParameter("ExpectNoise", "If true, expect noise hits", false);
	AddParameter("ConvolutionWidth", "If non-zero, convolve timing distribution with (approximately) a Gaussian of the given width", 0);
}

void
MMSLikelihood::Configure()
{
	GetParameter("InputPhotons", photons_name_);
	GetParameter("ExpectNoise", noise_);

	std::string spline_path;
	GetParameter("SplineTablePath", spline_path);
	spline_table_.read_fits(spline_path);

	//How does convolution work here?
	double convolution;
	GetParameter("ConvolutionWidth", convolution);
	if (convolution > 0) {
		// Use a triangular approximation to a Gaussian
		double knots[3] = {-2.*convolution, 0., 2.*convolution};
		spline_table_.convolve(0, knots, 3);
	}
}

void
MMSLikelihood::SetEvent(const I3Frame &fr)
{
	photons_ = fr.Get<I3PhotonSeriesMapConstPtr>(photons_name_);
	i3_assert(photons_);
	geo_ = fr.Get<I3GeometryConstPtr>(); // maybe assert
	i3_assert(geo_);
}

std::ofstream outfile("tracking.csv");
//outfile<<"pdf contribution"<< ","<<'dr'<<","<<'Ephi'<<","<<"dphi"<<","<<"tres"<<","<<"eX"<<","<<"eY"<<","<<"eZ"<<","<<"eT"<<","<<"eZenith"<<","<<"eAzimuth"<<'\n';


double
MMSLikelihood::GetLogLikelihood(const I3EventHypothesis &hypo)
{
	const I3Particle& part = *hypo.particle;
	
	
	double llh = 0;
	for(const auto& p : *photons_){ // for compressedphotonseries in compressedphotonseriesmap
		ModuleKey Mkey=p.first;
		OMKey om = OMKey(Mkey.GetString(), Mkey.GetOM(), 1);
		auto geo_it = geo_->omgeo.find(om);
		if(geo_it==geo_->omgeo.end())
			log_fatal_stream(om << " not found in geometry");
		
		
		
		
		I3Position cherenkov_emission_point;
		double t_direct=0, dist=0, Ephi=0, dphi=0;
		// cherenkovcalc only works for tracks
		// I3Calculator::CherenkovCalc(part,
		//     geo_it->second.position, cherenkov_emission_point,
		//     t_direct, dist, arrival_ang,
		//     1.35557, 1.34 /* fiducial group and phase n */);
		
		
		

		for(const auto& pulse : p.second){ // for photon in compressedphotonseries
			
			
			// calculate distance by hand
			I3Position diff =  geo_it->second.position - part.GetPos(); 

			dist += diff.Magnitude(); 
		
			t_direct += part.GetTime() + 1.34*dist / I3Constants::c;
			// Look at documentation to figure out syntax

			Ephi = acos((part.GetDir() * diff) / dist);
			// dataclasses.I3Direction Eangle = dataclasses.I3Direction([Ex, Ey, Ez])
			I3Direction photondir = pulse.GetDir();
			I3Position photonpos = pulse.GetPos();
			
			dphi = acos(std::max(-1.0, std::min(1.0,-1*(photondir*photonpos) / photonpos.Magnitude())));
			if (dphi < 0 or dphi > 1.5707963267948966){
				outfile<< "dphi: " << dphi<< '\n';
			}

			double splinecoords[4] = {dist, Ephi, dphi, pulse.GetTime() - t_direct}; 
			double pdf = spline_table_(splinecoords);
			
			if (!noise_) {
				if (pdf != 0) {
					llh += pdf; // Expects log pdf
					// outfile << '\n' << -pdf << ", "  << dist << ", " << Ephi << ", " << dphi << ", "<< pulse.GetTime() - t_direct << ", "<< part.GetPos().GetX() << ", " << part.GetPos().GetY() << ", " <<part.GetPos().GetZ() << ", " << part.GetTime() << ", " << part.GetDir().GetZenith() << ", " << part.GetDir().GetAzimuth();
					// std::cout<< pdf <<std::endl;
					continue;
				}
				else {
					llh += -30;
					// outfile << '\n' << 30 << ", "  << dist << ", " << Ephi << ", " <<dphi<<", " << pulse.GetTime() - t_direct << ", "<< part.GetPos().GetX() << ", " << part.GetPos().GetY() << ", " <<part.GetPos().GetZ() << ", " << part.GetTime() << ", " << part.GetDir().GetZenith() << ", " << part.GetDir().GetAzimuth();

				}
			}
			
		}
		//std::cout<< "llh: " << llh << std::endl;
	}
	//outfile << '\n' << -llh<< ' ' << part.GetPos().GetX() << ' ' << part.GetPos().GetY() << ' ' <<part.GetPos().GetZ() << ' ' << part.GetTime() << ' ' << part.GetDir().GetZenith() << ' ' << part.GetDir().GetAzimuth() << '\n';
	// Appears that Gulliver framework expects a negative number for llh, which it then inverts to minimize later. 
	return llh;
}

unsigned int
MMSLikelihood::GetMultiplicity()
{
	return 5000;
} // necessary to compile