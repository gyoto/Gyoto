gyoto=save(
           haveXerces=gyoto_haveXerces,
           haveCFITSIO=gyoto_haveCFITSIO,
           haveBoost=gyoto_haveBoost,
           haveUDUNITS=gyoto_haveUDUNITS,
           havePTHREAD=gyoto_havePTHREAD,
           haveMPI=gyoto_haveMPI,
           mpiInit=gyoto_mpiInit,
           mpiFinalize=gyoto_mpiFinalize,
           mpiInitialized=gyoto_mpiInitialized,
           mpiFinalized=gyoto_mpiFinalized,
           loadPlugin=gyoto_loadPlugin,

           Scenery=gyoto_Scenery,
           Scenery_rayTrace=gyoto_Scenery_rayTrace,
           Photon=gyoto_Photon,
           Metric=gyoto_Metric,
           Astrobj=gyoto_Astrobj,
           ThinDisk=gyoto_ThinDisk,
           Screen=gyoto_Screen,
           Spectrum=gyoto_Spectrum,

           debug=gyoto_debug,
           verbose=gyoto_verbose,

           Spectrometer=gyoto_Spectrometer,
           SpectroUniform=gyoto_SpectroUniform,
           SpectroComplex=gyoto_SpectroComplex,
           
           dontcatchSIGFPE=gyoto_dontcatchSIGFPE,
           dontcatchSIGSEGV=gyoto_dontcatchSIGSEGV,
           listRegister=gyoto_listRegister,

           is_Photon=is_gyoto_Photon,
           is_Astrobj=is_gyoto_Astrobj,
           is_Metric=is_gyoto_Metric,
           is_Spectrometer=is_gyoto_Spectrometer,
           is_Spectrum=is_gyoto_Spectrum,
           is_Screen=is_gyoto_Screen,
           is_Scenery=is_gyoto_Scenery,

           C=GYOTO_C,
           G=GYOTO_G,
           G_OVER_C_SQUARE=GYOTO_G_OVER_C_SQUARE,
           SUN_MASS=GYOTO_SUN_MASS,
           SUN_RADIUS=GYOTO_SUN_RADIUS,
           KPC=GYOTO_KPC
           
           );
  
func gyoto_namespace(args) {
  /* DOCUMENT gyoto_namespace, shortname=longname;
     If this version of yorick supports oxy object, the same as:
     save, gyoto, shortname=longname
     else, noop.
     SEE ALSO: oxy, save, restore
  */
  extern gyoto;
  keynames=args(-);
  for (i=1; i<=numberof(keynames); ++i) {
    save, gyoto, keynames(i), args(keynames(i)) ;
  }
}
wrap_args, gyoto_namespace;
