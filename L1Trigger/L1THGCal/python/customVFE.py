from L1Trigger.L1THGCal.hgcalVFEProducer_cfi import vfe_proc

def custom_old_digi(process):
    process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_ee.oldDigi = True
    process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_hesi.oldDigi = True
    process.hgcalVFEProducer.ProcessorParameters.linearizationCfg_hesc.oldDigi = True
    process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_ee.oldDigi = True
    process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_hesi.oldDigi = True
    process.hgcalVFEProducer.ProcessorParameters.calibrationCfg_hesc.oldDigi = True
    return process

def custom_hgcroc_oot(process,
                      oot_coefficients=vfe_proc.linearizationCfg_ee.oot_coefficients
                      ):
    parameters = vfe_proc.clone(
            linearizationCfg_ee = vfe_proc.linearizationCfg_ee.clone(oot_coefficients=oot_coefficients),
            linearizationCfg_hesi = vfe_proc.linearizationCfg_hesi.clone(oot_coefficients=oot_coefficients),
            linearizationCfg_hesc = vfe_proc.linearizationCfg_hesc.clone(oot_coefficients=oot_coefficients),
            )
    process.hgcalVFEProducer.ProcessorParameters = parameters
    return process


def custom_hgcroc_compression(process,
        exponentBits=vfe_proc.compressionCfg_ldm.exponentBits,
        mantissaBits=vfe_proc.compressionCfg_ldm.mantissaBits,
        rounding=vfe_proc.compressionCfg_ldm.rounding,
        truncationBits_ldm=vfe_proc.compressionCfg_ldm.truncationBits,
        truncationBits_hdm=vfe_proc.compressionCfg_hdm.truncationBits,
        ):
    parameters = vfe_proc.clone(
            compressionCfg_ldm = vfe_proc.compressionCfg_ldm.clone(
                exponentBits=exponentBits,
                mantissaBits=mantissaBits,
                truncationBits=truncationBits_ldm,
                rounding=rounding,
                ),
            compressionCfg_hdm = vfe_proc.compressionCfg_hdm.clone(
                exponentBits=exponentBits,
                mantissaBits=mantissaBits,
                truncationBits=truncationBits_hdm,
                rounding=rounding,
                ),
            )
    process.hgcalVFEProducer.ProcessorParameters = parameters
    return process
