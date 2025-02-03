# detector_geometry

Creates a csv file with CMS tracker module positions, rotations, thickness, bounds, pitches. 
An example csv is also in the repo. 
Most of the logic is defined in 

```trackerGeom/GeomDumper/plugins/GeomDumper.cc``` 

To run the project: 

```
cmsrel CMSSW_14_1_2
cd CMSSW_14_1_2/src
cmsenv
cmsRun trackerGeom/GeomDumper/test/geometryanalyzer_cfg.py
