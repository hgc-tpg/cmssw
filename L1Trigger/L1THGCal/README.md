# Documentation of the HGCAL TPG simulation
This README is meant to contain technical documentation on the HGCAL TPG simulation. 
Introductory, user-oriented documentation and installation recipes can be found in the HGCAL TPG Simulation [Twiki page](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HGCALTriggerPrimitivesSimulation).

## Code architecture
### Data Formats

### Producers, Processors and Implementations

### Standalone Emulators and Wrapping architecture

## Steps of the HGCAL TPG simulation
### Front-end: trigger path in the HGCROC

### Front-end: ECON-T

### Back-end: Clustering
#### Stage 1
#### Stage 2

### Back-end: Tower building


## Utilities
### Trigger geometry
The trigger geometry utility provides a unified interface to navigate between different geometrical objects and map hardware objects to the detector geometry. It comprises:
- Navigation between sensor cells and trigger cells
- Mapping between detector modules (hexaboards, tileboards) and trigger cells
- Mapping between modules and front-end optical links (lpGBT)
- Mapping between Stage 1 FPGAs and detector modules as well as lpGBTs
- Mapping between Stage 1 FPGAs and Stage 2 FPGAs

It also provides an interface to access the different underlying HGCAL geometries and topologies from the different sub-detectors (CE-E, CE-H-Si, CE-H-Sc), as well as some helper functions:
- Retrieve the number of elinks per detector module
- Retrieve the position of a trigger cell or of a detector module
- Retrieve active detector layers used in the trigger
- Check if a detector module is used in the trigger or disconnected (for instance if the corresponding detector layer is not used)

To be noted that previous versions of the trigger geometry utility implemented a way to retrieve trigger cell neighbors, which was used for 2D topological clustering. This feature is not available anymore in the currently existing versions.

Two versions of the trigger geometry utility are currently available:
- `HGCalTriggerGeometryV9Imp3` is the one used as default and implements all the link mappings between the front-end and the back-end, and within the back-end. The tag `V9` in the name means that it can be used with all detector geometries $\geq$ V9.
- `HGCalTriggerGeometryV9Imp2` is an older version in which the link mappings between the front-end and the back-end are not implemented.

They are stored in the [`geometries` plugin directory](plugins/geometries), and follow the interface defined in  [`interface/HGCalTriggerGeometryBase.h`](interface/HGCalTriggerGeometryBase.h). The trigger geometry utility is built by the ESProducer [`plugin/HGCalTriggerGeometryESProducer.cc`](plugins/HGCalTriggerGeometryESProducer.cc) so that it is made available through the `EventSetup`.

The parameters for the various trigger geometries are defined in [`python/customTriggerGeometry.py`](python/customTriggerGeometry.py), and the customization functions here can be used to switch from one trigger geometry to an other.
The default parameters of the trigger geometry are defined in [`python/l1tHGCalTriggerGeometryESProducer_cfi.py`](python/l1tHGCalTriggerGeometryESProducer_cfi.py) and alternative geometries can be used with the customization functions:
```python 
from L1Trigger.L1THGCal.customTriggerGeometry import custom_....
process = custom_....(process)
```

Link mappings accessed by the trigger geometry are provided as JSON files (`hgcal_trigger_link_mapping_*.json`) available in `L1Trigger/L1THGCal/data/`. Note that once integrated in CMSSW releases, these data files are stored on `cvmfs` under `$CMSSW_DATA_PATH/data-L1Trigger-L1THGCal` and can also be viewed from the Github [cms-data/L1Trigger-L1THGCal](https://github.com/cms-data/L1Trigger-L1THGCal) repository.


### Trigger tools
