#---------------------------------------------------------------------------
#
#  Makefile for executing the Donjon non-regression tests
#  Author : A. Hebert (2018-5-10)
#
#---------------------------------------------------------------------------
#
OS = $(shell uname -s | cut -d"_" -f1)
ifneq (,$(filter $(OS),SunOS AIX))
  MAKE = gmake
endif
ifeq ($(intel),1)
  fcompiler = intel
else
  fcompiler =
endif
all :
	$(MAKE) -C src
clean :
	$(MAKE) clean -C src
tests :
	./rdonjon 1Sap_vs_3Sap.x2m $(fcompiler) -quiet
	./rdonjon AFMtest.x2m $(fcompiler) -quiet
	./rdonjon Candu6.x2m $(fcompiler) -quiet
	./rdonjon candu6_iram.x2m $(fcompiler) -quiet
	./rdonjon Example1.x2m $(fcompiler) -quiet
	./rdonjon Example2.x2m $(fcompiler) -quiet
	./rdonjon FourChannels.x2m $(fcompiler) -quiet
	./rdonjon NCRtest.x2m $(fcompiler) -quiet
	./rdonjon SCR_RES_Interp.x2m $(fcompiler) -quiet
	./rdonjon SCWR64Core.x2m $(fcompiler) -quiet
	./rdonjon TC_EVO.x2m $(fcompiler) -quiet
	./rdonjon TC_XENON.x2m $(fcompiler) -quiet
	./rdonjon TINST.x2m $(fcompiler) -quiet
	./rdonjon Test4x4coupl_nompi.x2m $(fcompiler) -quiet
	./rdonjon Test_microdepl_NCR.x2m $(fcompiler) -quiet
	./rdonjon Test_IDET.x2m $(fcompiler) -quiet
	./rdonjon Test_PMAXS_Macro.x2m $(fcompiler) -quiet
	./rdonjon Test_PMAXS_Map.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH_void_diff_picard.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH_void.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH_void_bfgs.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH_void_lbfgs.x2m $(fcompiler) -quiet
	./rdonjon Test_SPH_void_newton.x2m $(fcompiler) -quiet
	./rdonjon Triga_Model3D.x2m $(fcompiler) -quiet
	./rdonjon VHTR_SPH.x2m $(fcompiler) -quiet
	./rdonjon VHTRc.x2m $(fcompiler) -quiet
	./rdonjon channel_mphy.x2m $(fcompiler) -quiet
	./rdonjon d2p_mco.x2m $(fcompiler) -quiet
	./rdonjon d2p_rfl.x2m $(fcompiler) -quiet
	./rdonjon d2p_sap.x2m $(fcompiler) -quiet
	./rdonjon fbr_core.x2m $(fcompiler) -quiet
	./rdonjon history_mphy.x2m $(fcompiler) -quiet
	./rdonjon pincell_mphy_thm.x2m $(fcompiler) -quiet
	./rdonjon pmox_sap.x2m $(fcompiler) -quiet
	./rdonjon pulseTHM_0d.x2m $(fcompiler) -quiet
	./rdonjon rep1300_reflector.x2m $(fcompiler) -quiet
	./rdonjon rep1300_reflector_sp1.x2m $(fcompiler) -quiet
	./rdonjon rep900.x2m $(fcompiler) -quiet
	./rdonjon rep900cluster.x2m $(fcompiler) -quiet
	./rdonjon rep900EnrichCOMPO.x2m $(fcompiler) -quiet
	./rdonjon rep900EnrichCOMPOhet.x2m $(fcompiler) -quiet
	./rdonjon rep900_iram.x2m $(fcompiler) -quiet
	./rdonjon rep900_mphy.x2m $(fcompiler) -quiet
	./rdonjon rep900_mphy_alt.x2m $(fcompiler) -quiet
	./rdonjon rep900_msap.x2m $(fcompiler) -quiet
	./rdonjon rep900_sim_defrod.x2m $(fcompiler) -quiet
	./rdonjon rep900_sim_recopy.x2m $(fcompiler) -quiet
	./rdonjon rep900_sim.x2m $(fcompiler) -quiet
	./rdonjon testNAPhet.x2m $(fcompiler) -quiet
	./rdonjon test_channBWR.x2m $(fcompiler) -quiet
	./rdonjon vver_test.x2m $(fcompiler) -quiet
