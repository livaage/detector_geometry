// -*- C++ -*-
//
// Package:    trackerGeom/GeomDumper
// Class:      GeomDumper
//
/**\class GeomDumper GeomDumper.cc trackerGeom/GeomDumper/plugins/GeomDumper.cc

 Description: Dumps the tracking module geometry into a csv file. 
 		Includes positions, rotations etc. 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Liv Helen Vage
//         Created:  Wed, 20 Nov 2024 13:37:15 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GeomDetType.h"
//#include "Geometry/TrackerGeometryBuilder/interface/GeomDetUnit.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometrySurface/interface/Bounds.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "FWCore/Utilities/interface/Exception.h" 
#include <iostream>
#include <fstream> 

//
// class declaration
//

class GeomDumper : public edm::stream::EDFilter<> {
public:
  explicit GeomDumper(const edm::ParameterSet&);
  ~GeomDumper() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  bool filter(edm::Event&, const edm::EventSetup&) override;
  void endStream() override; 
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
  const TrackerGeometry* geom_; 
  const edm::ESGetToken<TrackerTopology,TrackerTopologyRcd> topoToken_;
  const TrackerTopology* topo_;
  uint32_t find_volume(const TrackerTopology* topo, const DetId& id);
 // const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeometryToken_;
  //edm::Handle<TrackerDigiGeometryRecord> trackerGeom_; 
  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  //edm::EDGetTokenT<ExampleData> exampleToken_;
#endif
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GeomDumper::GeomDumper(const edm::ParameterSet& iConfig): 
  trackerGeometryToken_(esConsumes()), 
  topoToken_(esConsumes())
{
  //now do what ever initialization is needed
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  exampleToken_ = consumes<ExampleData>(iConfig.getParameter<edm::InputTag>("examples"));
#endif
#ifdef THIS_IS_EN_EVENTSETUP_EXAMPLE
  setupToken_ = esConsumes<SetupData, SetupRecord>();
#endif
}

GeomDumper::~GeomDumper() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

uint32_t GeomDumper::find_volume(const TrackerTopology* topo, const DetId& id){
	    if      ( topo == nullptr )                       return 0;
    else if ( id.subdetId()==2 && topo->side(id)==1 ) return 1; // IT endcap -ve
    else if ( id.subdetId()==1 && topo->side(id)==0 ) return 2; // IT barrel
    else if ( id.subdetId()==2 && topo->side(id)==2 ) return 3; // IT endcap +ve
    else if ( id.subdetId()==4 && topo->side(id)==1 ) return 4; // OT endcap -ve
    else if ( id.subdetId()==5 && topo->side(id)==0 ) return 5; // OT barrel
    else if ( id.subdetId()==4 && topo->side(id)==2 ) return 6; // OT endcap +ve
    else                                              return 0;
  };

// ------------ method called on each new Event  ------------
bool GeomDumper::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  std::ofstream outFile("cms_Extended2026D110_ModulePositions.csv"); 
  if (!outFile) {
    std::cerr << "Error opening file for writing." << std::endl;
    }
  outFile << "volume_id,layer_id,module_id,raw_id,cx,cy,cz,rot_xu,rot_xv,rot_xw,rot_yu,rot_yv,rot_yw,rot_zu,rot_zv,rot_zw,module_t,module_minhu,module_maxhu,module_hv,pitch_u,pitch_v\n"; 
#ifdef THIS_IS_AN_EVENT_EXAMPLE
  //ExampleData const& in = iEvent.get(exampleToken_);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //SetupData const& setup = iSetup.getData(setupToken_);
#endif
//      edm::ESHandle<TrackerGeometry> trackerGeometry;
//    iSetup.get<TrackerDigiGeometryRecord>().get(trackerGeometry);
	geom_ = &iSetup.getData(trackerGeometryToken_);   
	topo_ = &iSetup.getData(topoToken_);

//const TrackerGeometry& trackerGeometry = iSetup.getData(trackerGeometryToken_);
for (const auto& det : geom_->dets()) {
        const GeomDet* geomDet = det;
	const TrackerTopology* topo = topo_; 
        DetId detId = geomDet->geographicalId();
	int subdetId = detId.subdetId(); 
	uint32_t volume = find_volume(topo_, detId); 
        // Volume, layer, module IDs
 //       int volume_id = detId.subdetId();    // Example: Subdetector ID
 //       int layer_id = detId.layer();       // Assuming layer ID is retrievable
        int module_id = detId.rawId();      // Unique module ID
	int mod_id = topo->module(detId);
	int layer_id = topo->layer(detId); 
	//int layer = (module_id >> 3) & 0x7;
	if (subdetId == PixelSubdetector::PixelBarrel){
		int module_id = (module_id >> 9) & 0xF; 
	}else{
		int module_id = 0; 
	}
	

        // Position of the local origin in global coordinates
        const GlobalPoint& globalPosition = geomDet->position();
        double cx = globalPosition.x();
        double cy = globalPosition.y();
        double cz = globalPosition.z();

        // Rotation matrix from local to global coordinates
        const Surface::RotationType& rotation = geomDet->surface().rotation();
        double rot_xu = rotation.xx();
        double rot_xv = rotation.xy();
        double rot_xw = rotation.xz();
        double rot_yu = rotation.yx();
        double rot_yv = rotation.yy();
        double rot_yw = rotation.yz();
        double rot_zu = rotation.zx();
        double rot_zv = rotation.zy();
        double rot_zw = rotation.zz();

		// Module dimensions
        const Bounds* bounds = &(geomDet->surface().bounds());
	const RectangularPlaneBounds* rectBounds = dynamic_cast<const RectangularPlaneBounds*>(bounds);
	double module_minhu = 0; 
	double module_maxhu = 0; 
	double module_hv = 0; 
	if (rectBounds) {
    	module_minhu = rectBounds->width() / 2.0;  // Minimum half-length (local u)
    	module_maxhu = rectBounds->width() / 2.0;   // Maximum half-length (local u)
	module_hv = rectBounds->length() / 2.0;  // Half-length in local v direction
	}
        //double module_minhu = bounds->width() / 2.0;
        //double module_maxhu = bounds->width() / 2.0;  // Assuming uniform bounds
        //double module_hv = bounds->length() / 2.0;
        double module_t = rectBounds->thickness() / 2.0;
	double pitch_u = 0; 
	double pitch_v = 0; 

if (subdetId < 4) {
    //std::cout << "DetUnit type: " << typeid(*detUnit).name() << std::endl;
    const GeomDetUnit* detUnit = geom_->idToDetUnit(detId);
//if (detId.subdetId() == 1 || detId.subdetId() == 2) { // Pixel Barrel or Pixel Endcap
   std::cout<<"in det"<<std::endl;
    const PixelGeomDetUnit* pixelDet = dynamic_cast<const PixelGeomDetUnit*>(detUnit);
    if (pixelDet) {
        std::cout << "This is a pixel module!" << std::endl;
        const PixelTopology& topology = pixelDet->specificTopology();
        std::pair<float, float> pitchUV = topology.pitch(); // Pitch in (u, v)
        pitch_u = pitchUV.first;  // Pitch in u direction
        pitch_v = pitchUV.second; // Pitch in v direction
        std::cout << "Pitch (u, v): (" << pitch_u << ", " << pitch_v << ")" << std::endl;
    } else {
        std::cout << "Not a pixel module." << std::endl;
    }
} //else {
   // std::cout << "Invalid DetUnit for detId: " << detId.rawId() << std::endl;
  //  continue;
//}

//}

        // Pitch sizes (can be retrieved from GeomDetType or associated data)

	outFile << volume << "," << layer_id << "," << mod_id << "," << module_id <<"," << cx << "," <<  cy << "," << cz << "," << rot_xu << "," << rot_xv << "," << rot_xw << "," << rot_yu 
		<< "," << rot_yv << "," << rot_yw << "," << rot_zu << "," << rot_zv <<"," << rot_zw << "," << module_t <<"," <<  module_minhu << "," << module_maxhu << "," << module_hv << "," << pitch_u << "," << pitch_v << "\n"; 

	//For debugging 

        // Print or store the information
        std::cout << ", Module ID: " << module_id
                  << ", Position: (" << cx << ", " << cy << ", " << cz << ")"
                  << ", Rotation: [" << rot_xu << ", " << rot_xv << ", " << rot_xw
                  << " | " << rot_yu << ", " << rot_yv << ", " << rot_yw
                  << " | " << rot_zu << ", " << rot_zv << ", " << rot_zw << "]"
                  << ", Thickness: " << module_t
                  << ", Bounds: [" << module_minhu << ", " << module_maxhu << ", " << module_hv << "]"
                  << ", Pitches: [" << pitch_u << ", " << pitch_v << "]"
		  <<", Volume: "<<volume<<" , "
		  <<", Layer: "<<layer_id<<" , "
		  <<", Module: "<<module_id<<" , "
                  << std::endl;
    }
    outFile.close(); 
//      for (const auto& det : geom_->dets()) {
//        if (det) {
//		std::cout<<"Detector position "<<det->position()<<std::endl;
        //    edm::LogInfo("GeomDumper") << "Detector position: " << det->position();
 //      }
  //  }

  return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void GeomDumper::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void GeomDumper::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
GeomDumper::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
GeomDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GeomDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GeomDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GeomDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GeomDumper);
