# Stage-2 semi-emulator

This directory contains the older hybrid Stage-2 implementation used by
`HGCalBackendLayer2ProcessorSemiEmulator`.  It performs floating-point
triangular clustering and fixed-width cluster-property calculations.  It is a
validation reference and does not correspond to a released firmware design.

The code was imported from `indra-ehep/hgcal-tpg-fe` commit
`6c1806002cf47f278a6a34f9e7a7096f64500cf5`, through the integration in
`EmyrClement/cmssw` commit `c1a9cc5b3ec482c461c07e7a0d2be6c807ea1e60`.
The packed cluster type is kept local to this implementation instead of being
added to the persistent CMSSW data model.  `TPGTCFloats` was also changed to
use call-local rotation temporaries so concurrent CMSSW streams do not share
mutable state.
