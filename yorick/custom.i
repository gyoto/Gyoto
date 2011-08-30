plug_dir,_(["./", "./stdplug/"], plug_dir());
set_path, "stdplug:"+get_path();
//#include "check.i"

/* This should be the final line in your custom.i file-- it implements
   the default command line processing (see help for process_argv).  */
command_line= process_argv();
