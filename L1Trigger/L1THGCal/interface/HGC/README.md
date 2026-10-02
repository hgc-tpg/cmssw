# Imported Stage-2 emulator core

This directory is a mechanically imported snapshot of `hgcal-tpg/Stage2`,
commit `1a22402`, from `Emulator/HGC`.  Its source of truth remains the Stage2
repository until the common-core migration is complete.

The adapter in `plugins/backend` supplies trigger-cell geometry and the
temporary virtual transport mapping.  The import has a small set of
intentional integration changes: CMSSW-qualified includes, heap-backed
workspace storage for the cluster pipeline, safe handling of empty decoder
inputs, a single-cell decoder entry point, and corrections to the fixed-point
scales used by the example logical chain.  The latter corrections remain to
be confirmed against the firmware developer's full-chain conventions. Keep
future CMSSW-specific logic in the adapter.

To update this snapshot, begin with a mechanical import of the same Stage2
subdirectory, then reapply those small compatibility changes and record the
new source commit above.  The longer-term preferred arrangement is to extract
this core to a versioned common package, consumed by both the firmware and
CMSSW wrappers.
