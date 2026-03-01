cet_make_library(LIBRARY_NAME OscLib
  SOURCE
    EarthModel.cxx
    IOscCalc.cxx
    IOscCalcSterile.cxx
    OscCalc.cxx
    OscCalcAnalytic.cxx
    OscCalcCPT.cxx
    OscCalcDMP.cxx
    OscCalcDumb.cxx
    OscCalcGeneral.cxx
    OscCalcNuFast.cxx
    OscCalcPMNS_CPT.cxx
    OscCalcPMNS_NSI.cxx
    OscCalcPMNS.cxx
    OscCalcPMNSOpt.cxx
    OscCalcPMNSOptEigen.cxx
    OscCalcSterile.cxx
    OscCalcSterileBeam.cxx
    OscCalcSterileEigen.cxx
    PMNS_DMP.cxx
    PMNS_NSI.cxx
    PMNS_Sterile.cxx
    PMNS.cxx
    PMNSOpt.cxx

  LIBRARIES
    OscLib_options
)

install_headers()
install_source()
