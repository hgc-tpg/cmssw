# Documentation of the HGCAL TPG simulation
This README is meant to contain technical documentation on the HGCAL TPG simulation. 
Introductory, user-oriented documentation and installation recipes can be found in the HGCAL TPG Simulation [Twiki page](https://twiki.cern.ch/twiki/bin/viewauth/CMS/HGCALTriggerPrimitivesSimulation).

## Code architecture
### Data Formats

### Producers, Processors and Implementations

### Standalone Emulators and Wrapping architecture

## Steps of the HGCAL TPG simulation
### Front-end: trigger path in the HGCROC
The simulation of the trigger path in the HGCROC (called also `VFE` in the code, for very-front-end) performs the following tasks:
- Linearisation of the input digitized charges (putting the ADC and TOT values to the same linear scale), taken from the HGCAL `digis`
- Forming trigger cells (sums of 4 or 9 sensor cells) from linearized input charges
- Compression of the trigger cell energies on a floating point format

These steps are called from [`plugins/veryfrontend/HGCalVFEProcessorSums.cc`](plugins/veryfrontend/HGCalVFEProcessorSums.cc), and configured with [`python/l1tHGCalVFEProducer_cfi.py`](python/l1tHGCalVFEProducer_cfi.py). This configuration can be customized with customization functions available in [`python/customVFE.py`](python/customVFE.py). The actual implementations of the processing steps are stored in [`src/veryfrontend`](src/veryfrontend).

:warning: To be noted that at the moment the calibration of the trigger cell charges to energies in GeV is called from there while in reality it will be done in the ECON-T.

### Front-end: ECON-T
The ECON-T simulation implements several data reduction strategies applied on the trigger cells sent from the HGCROC trigger path. The data reduction strategies available are the following:
- `Threshold`: Selecting a variable subset of trigger cells passing a threshold on their transverse energy. This is the one used as default.
- `BestChoice`, `BC`: Selecting a fixed number of trigger cells per detector module. The trigger cells with the highest transverse energies are retained. 
- `SuperTriggerCell`, `STC`: Aggregating further trigger cells into coarser objects called Super Trigger Cells. Several flavours of Super Trigger Cells are available and configurable.
- `AutoEncoder`: Applying trigger cell data compression with an auto-encoder neural network.

Different data reduction strategies can be configured in the different sub-detectors (CE-E, CE-H-Si, CE-H-Sc), and in particular `BestChoice` can be used in the CE-E with `SuperTriggerCell` in the CE-H (which corresponds to the so called `BC+STC` algorithm).

The ECON-T simulation also builds module sums (called `trigger sums` in the code), which are energy sums over entire detector modules. All trigger cells, or only those not selected by the selection algorithms described above, can be summed in module sums.

These steps are called from [`plugins/concentrator/HGCalConcentratorProcessorSelection.cc`](plugins/concentrator/HGCalConcentratorProcessorSelection.cc), and configured with [`python/l1tHGCalConcentratorProducer_cfi.py`](python/l1tHGCalConcentratorProducer_cfi.py). This configuration can be customized with customization functions available in [`python/customTriggerCellSelect.py`](python/customTriggerCellSelect.py) and [`python/customTriggerSums.py`](python/customTriggerSums.py):
```python
from L1Trigger.L1THGCal.customTriggerCellSelect import custom_....
from L1Trigger.L1THGCal.customTriggerSums import custom_....
process = custom_....(process)
```
The actual implementations of the algorithms are stored in [`src/concentrator`](src/concentrator).

:warning: To be noted that the calibration step from charges to energies in GeV is currently called from the HGCROC trigger path simulation, while it is in reality done in the ECON-T.

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
