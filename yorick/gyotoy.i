#!/usr/bin/env yorick -i
/*
    Copyright 2007 F. Rigaut
    Copyright 2011, 2013 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
   gyotoy.i
   main function to call the pygtk GUI to gyotoy.
   syntax: yorick -i gyotoy.i
*/

require,"pyk.i";
require, "gyoto.i";
require, "gyoto_std.i";
require, "pl3d.i";

extern GYOTOY_PYTHON3_DEFAULT;
extern GYOTOY_PYTHON3;
/* DOCUMENT extern GYOTOY_PYTHON3_DEFAULT;

     Gyotoy tries hard to find a working python3 installation.

     Users can set GYOTOY_PYTHON3 to there preferred version of
     python3, using either the simple command name or full path. This
     version of python3 must have PyGObject installed, and the GObject
     Introspection files for Gtk3 must be installed on the system.

     GYOTOY_PYTHON3 can be set in the environment, for instance in
     .profile:
      export GYOTOY_PYTHON3=/usr/local/bin/python3.very.experimental
     or in Yorick, for instance in a file under Y_USER/i-start/:
      GYOTOY_PYTHON3 = "/usr/local/bin/python3.very.experimental";
      
     Note that Y_USER represents the user's main Yorick
     directory. Type Y_USER at the Yorick prompt to check it. It's
     usually one of ~/.yorick, ~/Yorick and ~/yorick. See
     Y_SITE/i/custom.i for further details.
     
     Distributors may set GYOTOY_PYTHON3_DEFAULT to the Right value
     for their system, for instance /usr/bin/python3.2mu or
     /opt/local/bin/python3.3. Users will by default use this
     executable while still being able to customize it using the
     GYOTOY_PYTHON3 variable. GYOTOY_PYTHON3_DEFAULT can be set only
     in Yorick, typically by dropping a file in Y_SITE/i-start.

   SEE ALSO: Y_SITE

 */

extern _gyotoy_running;
extern _gyotoy_reticle;
extern _gyotoy_wid;
extern _gyotoy_wstyle;
extern _gyotoy_filename;
extern _gyotoy_stand_alone;
extern _gyotoy_particle_is_massless;
extern _gyotoy_metric;
extern _gyotoy_particle;
extern _gyotoy_t1;
extern _gyotoy_initcoord;
extern _gyotoy_mass;
extern _gyotoy_distance;
extern _gyotoy_unit;
extern _gyotoy_txyz;
extern _gyotoy_particle_to_load;
extern _gyotoy_parent_id;
extern _gyotoy_inhibit_redraw;
extern _gyotoy_metric_file;
extern _gyotoy_nsteps;
/* DOCUMENT extern _gyotoy_*;
     Some variable holding state information for use inside gyotoy.
     _running: true if GTK interface is running;
     _reticle: true if reticle should be displayed;
     _wid: Yorick window ID of the Gyotoy plot area;
     _wstyle: window style (normally "nobox.gs");
     _gyotoy_filename: name of the last file read (careful: can be
     the py file itself)
     _gyotoy_stand_alone: true if the user is not supposed to ever see
     the Yorick prompt.
     _gyotoy_particle_is_massless: true for a photon.
     _gyoto_metric: the actual metric object.
   SEE ALSO: gyotoy
 */
if (is_void(_gyotoy_nsteps)) _gyotoy_nsteps=100 ;

func gyotoy_reset(void){
  extern _gyotoy_running;
  extern _gyotoy_reticle;
  extern _gyotoy_wid;
  extern _gyotoy_wstyle;
  extern _gyotoy_filename;
  extern _gyotoy_stand_alone;
  extern _gyotoy_particle_is_massless;
  extern _gyotoy_metric;
  extern _gyotoy_particle;
  extern _gyotoy_t0;
  extern _gyotoy_t1;
  extern _gyotoy_initcoord;
  extern _gyotoy_mass;
  extern _gyotoy_distance;
  extern _gyotoy_unit;
  extern _gyotoy_txyz;
  extern _gyotoy_particle_to_load;

  _gyotoy_particle_is_massless=0;
  _gyotoy_metric=gyoto_KerrBL(spin=0.995);
  _gyotoy_initcoord=[0.,10.791,pi/2,0.,0.,0.,0.016664];
  _gyotoy_t1=3000.;
  _gyotoy_t0=0.;
  _gyotoy_particle=gyoto_Star(metric=_gyotoy_metric,
                              initcoord=_gyotoy_initcoord(1:4),
                              _gyotoy_initcoord(5:7),
                              xfill=_gyotoy_t1);
  _gyotoy_txyz=_gyotoy_particle(get_txyz=1);
  _gyotoy_mass=4e6;
  _gyotoy_distance=8.;
  _gyotoy_unit="geometrical";
}
gyotoy_reset;

func gyotoy_set_unit(unit) {
  extern _gyotoy_unit;
  _gyotoy_unit=unit;
  gyotoy_redraw;
}

func gyotoy_set_t1(t1) {
  extern _gyotoy_txyz, _gyotoy_t1;
  /*  if ( (_gyotoy_t1 > _gyotoy_t0 && t1<_gyotoy_t1) ||
       (_gyotoy_t1 < _gyotoy_t0 && t1>_gyotoy_t1) )
       _gyotoy_particle, reset=;*/
  _gyotoy_t1=t1;
  if (catch(0x08)) return; // avoid breaking in case of v>c
  //  _gyotoy_txyz=_gyotoy_particle(xfill=t1, get_txyz=1);
  gyotoy_compute_and_draw;
}

func gyotoy_set_distance(dist) {
  extern _gyotoy_distance;
  _gyotoy_distance=dist;
  gyotoy_redraw;
}

func gyotoy_set_KerrBL_metric( spin ) {
  extern _gyotoy_metric, _gyotoy_particle, _gyotoy_txyz;
  _gyotoy_metric = gyoto_KerrBL ( spin = spin );
  if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_particle,metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_set_metric( fname ) {
  extern _gyotoy_metric, _gyotoy_particle, _gyotoy_txyz, _gyotoy_metric_file;
  if (catch(0x08)) {
    // avoid breaking in case bad file
    gyotoy_warning, "Unable to load metric. Is this a GYOTO XML description file?";
    return;
  }
  if (is_string(fname)) metric=gyoto_Metric(fname);
  if (typeof(metric)!="gyoto_Metric") gyotoy_warning, "Failed to set metric";
  _gyotoy_metric=metric;
  _gyotoy_metric_file=fname;
  if (catch(0x08)) {
    // avoid breaking in case of v>c or other problem
    gyotoy_warning, "metric loaded but orbit computation failed";
    return;
  }
  _gyotoy_particle,metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_set_initcoord(t0, r0, theta0, phi0,
                          rprime0, thetaprime0, phiprime0) {
  extern _gyotoy_initcoord, _gyotoy_txyz, _gyotoy_inhibit_redraw;
  _gyotoy_initcoord=[t0, r0, theta0, phi0,
                     rprime0, thetaprime0, phiprime0];
  _gyotoy_t0=t0;
  if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_particle, initcoord=_gyotoy_initcoord(1:4),_gyotoy_initcoord(5:7);
  if (_gyotoy_inhibit_redraw) return;
  gyotoy_compute_and_draw;
}

func gyotoy_redraw(void){
/* DOCUMENT gyotoy_redraw

    Redraws the Gyotoy plot:
     - convert to UNIT;
     - erase the plot;
     - if requested, draws the reticle;
     - plot the orbit (which triggers drawing the 3d cage).
     
   SEE ALSO: gyoto_Kerr_orbit, gyoto_convert, gyoto_reticle,
             gyoto_plg3
 */
  //  [a, e, l, q];
  //  delta; long(niter);
  // Purpose of the "if": gyoto_Kerr_orbit will return 0 if v>c

  extern _gyotoy_txyz, _gyotoy_inhibit_redraw;
  if (_gyotoy_inhibit_redraw) return;
  if (numberof(_gyotoy_txyz)>1) {
    t=_gyotoy_txyz(,1);
    n0=abs(t-_gyotoy_t0)(mnx);
    n1=abs(t-_gyotoy_t1)(mnx);
    if (n1>n0) ++n1; else --n1;
    if (n1<1) n1=1;
    if (n1>numberof(t)) n1=numberof(t);
    if (n1<n0) {
      n2=n0;
      n0=n1;
      n1=n2;
    }
    x=_gyotoy_txyz(n0:n1,2);
    y=_gyotoy_txyz(n0:n1,3);
    z=_gyotoy_txyz(n0:n1,4);
    gyoto_convert, x, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    gyoto_convert, y, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    gyoto_convert, z, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    rmax=sqrt(max(x^2+y^2+z^2));

    window, _gyotoy_wid;
    clear3;
    if (_gyotoy_reticle) gyoto_reticle, rmax*1.41421;
    if (gyotoy_origin) gyoto_plmksky,0,0;
    gyoto_plg3, x, y, z;
    limit3,-rmax,rmax,-rmax,rmax,-rmax,rmax;
  }

}

func gyotoy_toggle_reticle(void) {
/* DOCUMENT gyotoy_toggle_reticle
    Toggles sky reticle on/off
   SEE ALSO: gyoto_reticle
 */
  extern _gyotoy_reticle;
  if (is_void(_gyotoy_reticle)) _gyotoy_reticle=1;
  else _gyotoy_reticle=1-_gyotoy_reticle;
}

func gyotoy_window_init(parent_id)
// initialize window in GTK frontend
{
  extern _gyotoy_wid, _gyotoy_wstyle, _gyotoy_filename, _gyotoy_particle, _gyotoy_parent_id, _gyotoy_inhibit_redraw;
  if (_gyotoy_parent_id==parent_id) return;
  _gyotoy_parent_id=parent_id;
  ok=pyk("sleep(0.1)")
  _gyotoy_wstyle="nobox.gs";
  window,_gyotoy_wid,wait=1,parent=parent_id,style="nobox.gs",dpi=90;
  limits, square=1;
  gnomon,1;
  cage3,1;
  //orient3,pi,-pi/2;
  pldefault, marks=0;
  if (!is_void(_gyotoy_particle_to_load)) gyotoy_set_particle,_gyotoy_particle_to_load;
  else if (_gyotoy_filename) gyotoy_import,_gyotoy_filename;
  else pyk,"compute_orbit('rien')";
  pyk,"orient3('rien')";
  pyk,"builder.get_object('metric_type').set_active(0)";
  _gyotoy_inhibit_redraw=0;
  gyotoy_compute_and_draw;
}

func gyotoy_toggle_window_style(parent_id) {
// toggle between nobox.gs and work.gs. Currently unused.
  extern _gyotoy_wid, _gyotoy_wstyle;
  winkill,_gyotoy_wid;
  if (_gyotoy_wstyle=="nobox.gs") _gyotoy_wstyle="work.gs";
  else _gyotoy_wstyle="nobox.gs";
  _gyotoy_wstyle;
  window,_gyotoy_wid,wait=1,parent=parent_id,style=_gyotoy_wstyle;
}

func gyotoy_quit(void) {
// called when GTK window is closed
  extern _gyotoy_running;
  _gyotoy_running=0;
  _gyotoy_parent_id=0;
  if (_gyotoy_stand_alone) quit;
}

func gyotoy_checkvers_cb(msg) {
  extern _gyotoy_pyvers;
  _gyotoy_pyvers=msg;
}

func gyotoy(filename) {
/* DOCUMENT gyotoy [,filename]
         or gyotoy, star
         or gyotoy, photon
     Launch Gyotoy GTK interface
   SEE ALSO:
 */
  extern _pyk_proc, _gyotoy_filename, _gyotoy_running, _gyotoy_particle;
  extern _gyotoy_particle_to_load, _gyotoy_parent_id, _gyotoy_inhibit_redraw;
  _gyotoy_running=1;
  _gyotoy_filename=[];
  _gyotoy_particle_to_load=[];
  _gyotoy_parent_id=0;
  _gyotoy_inhibit_redraw=1;
  
  python_script=find_in_path("gyotoy.py", takefirst=1,
                             path=pathform(_(get_cwd(),
                                             _(Y_SITES,
                                               Y_SITE)+"python/")));
 
  if (strpart(python_script, 1:2)=="~/")
    python_script=get_home()+"/"+strpart(python_exec, 3:);
  
  glade_file= find_in_path("gyotoy.xml", takefirst=1,
                           path=pathform(_(get_cwd(), _(Y_SITES,Y_SITE)+"glade/")));

  if (strpart(glade_file, 1:2)=="~/")
    glade_file=get_home()+"/"+strpart(glade_file, 3:);

  gyotoytop = dirname(glade_file)+"/";

  //// SETTING GYOTOY_PYTHON3
  // First check if it has been set by the user
  if (is_void(GYOTOY_PYTHON3)) GYOTOY_PYTHON3 = get_env("GYOTOY_PYTHON3");
  // Then use distribution default
  if (!GYOTOY_PYTHON3) GYOTOY_PYTHON3=GYOTOY_PYTHON3_DEFAULT;
  // Then look for standard names
  PATH=pathsplit(get_env("PATH"));
  for (p=1; p<=numberof(PATH); ++p) {
    if (strpart(PATH(p), 0:0) != "/") PATH(p)+="/";
  }
  PATH=pathform(PATH);
  // Try python3
  if (is_void(GYOTOY_PYTHON3))
    GYOTOY_PYTHON3 = find_in_path("python3", takefirst=1, path=PATH);

  // Try python, check whether it's version 3
  if (is_void (GYOTOY_PYTHON3)) {
    GYOTOY_PYTHON3 = find_in_path("python", takefirst=1, path=PATH);
    if (!is_void(GYOTOY_PYTHON3)) {
        extern _gyotoy_pyvers;
        _gyotoy_pyvers="";
        pyproc=spawn([GYOTOY_PYTHON3, "-V"], noop ,gyotoy_checkvers_cb);
        count=0;
        while (_gyotoy_pyvers=="" && count < 50) {
          pause, 100;
          ++count;
        }
        if (strpart(_gyotoy_pyvers, 1:8) != "Python 3")
          GYOTOY_PYTHON3 = [];
    }
  }
  
  // Check python3.m for some values of m
  for (minor=10; minor>=0 && is_void(GYOTOY_PYTHON3); --minor) {
    GYOTOY_PYTHON3 = find_in_path("python3."+pr1(minor),
                                  takefirst=1, path=PATH);
  }
  if (is_void(GYOTOY_PYTHON3)) error, "Cannot find python3 executable";

  // build command to spawn:
  pyk_cmd=[GYOTOY_PYTHON3, python_script, gyotoytop];
  
  // spawn it and attach to _pyk_callback (see pyk.i):
  // that starts python and pass it the path to the glade file
  if (is_string(filename)) _gyotoy_filename=filename;
  else if ((is_gyoto_Astrobj(filename) && filename(kind=)=="Star")
           ||typeof(filename)=="gyoto_Photon")
    _gyotoy_particle_to_load=filename;
  else _gyotoy_particle_to_load=_gyotoy_particle;
  _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  if (is_void(_pyk_proc)) error, "Failed to launch python";
  if (_gyotoy_filename &&
      (strglob("*.dat",_gyotoy_filename,case=0) ||
       strglob("*.txt",_gyotoy_filename,case=0) ||
       strglob("*.xml",_gyotoy_filename,case=0) ))
    pyk,"set_filename('"+_gyotoy_filename+"')";
}

func gyotoy_set_particle(part) {
  extern _gyotoy_particle, _gyotoy_metric, _gyotoy_initcoord, _gyotoy_txyz,
    _gyotoy_metric_file, _gyotoy_filename, _gyotoy_inhibit_redraw;

  rdr=_gyotoy_inhibit_redraw;
  _gyotoy_inhibit_redraw=1;
  
  omtype=_gyotoy_metric(kind=);
  oldmass=_gyotoy_particle_is_massless;
  
  _gyotoy_particle=part;
  _gyotoy_metric = part(metric=);
  _gyotoy_initcoord=part(initcoord=);
  
  if (is_gyoto_Astrobj(part)) {
    part_type="star";
    _gyotoy_particle_is_massless=0;
  }
  else {
    part_type="photon";
    _gyotoy_particle_is_massless=1;
  }
  if (oldmass != _gyotoy_particle_is_massless)
    ok=pyk("set_parameter('particle_type','"+part_type+"')");

  // Metric & projection
  if ((mtype=_gyotoy_metric(kind=))=="KerrBL") {
    if (omtype != "KerrBL")
      ok=pyk("set_parameter('metric_type', 'kerrbl')");
    ok=pyk("set_parameter('spin',"+
           swrite(format="%.12f",_gyotoy_metric(spin=))+")");
  } else {
    if (otype == "KerrBL")
      ok=pyk("set_parameter('metric_type', 'file')");
    if (_gyotoy_filename)
      ok=pyk("set_parameter('metric_file', '"+_gyotoy_filename+"')");
    _gyotoy_metric_file=_gyotoy_filename;
  }
  
  //ok=pyk("set_parameter('incl',"+swrite(format="%.12f",metric(get_inclination=1)*rad2deg)+")");
  //ok=pyk("set_parameter('paln',"+swrite(format="%.12f",metric(get_paln=1)*rad2deg)+")");
  //ok=pyk("set_parameter('phase',"+swrite(format="%.12f",metric(get_argument=1)*rad2deg)+")");
  //ok=pyk("set_parameter('distance',"+swrite(format="%.12f",metric(get_distance=1))+")");
  //  m_sun = 1.98843e30;     // kg
  //  ok=pyk("set_parameter('mass',"+swrite(format="%.12f",_gyotoy_metric(mass=)/m_sun)+")");

  // Initial condition
  coord = part(initcoord=);
  ok=pyk("set_parameter('t0',"+swrite(format="%.12f",coord(1))+")");
  ok=pyk("set_parameter('r0',"+swrite(format="%.12f",coord(2))+")");
  ok=pyk("set_parameter('theta0',"+swrite(format="%.12f",coord(3))+")");
  ok=pyk("set_parameter('phi0',"+swrite(format="%.12f",coord(4))+")");
  if (part_type=="star") fact=1./coord(5);
  else fact=1.
  ok=pyk("set_parameter('rprime0',"+swrite(format="%.12f",coord(6)*fact)+")");
  ok=pyk("set_parameter('thetaprime0',"+swrite(format="%.12f",coord(7)*fact)+")");
  ok=pyk("set_parameter('phiprime0',"+swrite(format="%.12f",coord(8)*fact)+")");

  // Wait for parameters to have reached glade
  //    pyk,"compute_orbit('rien')";

  // bug ?
  if (mtype != "KerrBL")
    ok=pyk("set_parameter('metric_type', 'file')");

  _gyotoy_inhibit_redraw=rdr;
  gyotoy_compute_and_draw;
  
}

func gyotoy_set_nsteps(nsteps) {
  extern _gyotoy_nsteps;
  nsteps=long(nsteps);
  if (nsteps <=0) nsteps=1;
  _gyotoy_nsteps=nsteps;
}

func gyotoy_compute_and_draw(rien) {
  
  extern _gyotoy_particle, _gyotoy_redrawing, _gyotoy_cancel, _gyotoy_nsteps;
  extern _gyotoy_t1, _gyotoy_inhibit_redraw;
  if (_gyotoy_inhibit_redraw) return;
  
  if (_gyotoy_redrawing) {
    _gyotoy_cancel=1;
    return;
  }
  _gyotoy_cancel=0;

  if (is_void(_gyotoy_redrawing)) _gyotoy_redrawing=0;
  ++_gyotoy_redrawing;

  t0 = _gyotoy_particle(initcoord=)(1);
  _gyotoy_txyz=_gyotoy_particle(get_txyz=1);
  if (_gyotoy_t1>=t0)
    t  = _gyotoy_txyz(0,1);
  else
    t  = _gyotoy_txyz(1,1);
  dt = (_gyotoy_t1-t0)/_gyotoy_nsteps;
  if (dt==0) {
    gyotoy_warning, "t0 and t1 too close.";
    --_gyotoy_redrawing;
    return;
  }
  if ( ((t-t0)/(_gyotoy_t1-t0)) >=1 ){
    pause, 10;
    pyk,"set_play_image('gtk-media-play')";
    pyk,"set_fraction("+pr1((t-t0)/(_gyotoy_t1-t0))+")";
    gyotoy_redraw;
    --_gyotoy_redrawing;
    return;
  }

  pyk,"set_fraction("+pr1((t-t0)/(_gyotoy_t1-t0))+")";
  pyk,"set_play_image('gtk-media-pause')";

  dir = sign(_gyotoy_t1-t0);
  
  for ( ;
        ((dir == 1 && t<=_gyotoy_t1+dt) ||
         (dir ==-1 && t>=_gyotoy_t1+dt) )
          && !_gyotoy_cancel;
        t+=dt) {
    _gyotoy_txyz=_gyotoy_particle(xfill=t, get_txyz=1);
    if ( ( (dir== 1) && _gyotoy_txyz(0,1) < t) ||
         ( (dir==-1) && _gyotoy_txyz(1,1) > t) )
      t = _gyotoy_t1+2*dt;
    gyotoy_redraw;
    pause,10;
    pyk,"set_fraction("+pr1((t-t0)/(_gyotoy_t1-t0))+")";
  }
  pause, 10;
  pyk,"set_play_image('gtk-media-play')";
  pyk,"set_fraction("+pr1((t-t0)/(_gyotoy_t1-t0))+")";
  gyotoy_redraw;
  --_gyotoy_redrawing;
}

func gyotoy_inhibit_redraw(mode) {
  extern _gyotoy_inhibit_redraw;
  _gyotoy_inhibit_redraw=mode;
}

func gyotoy_export(filename) {
/* DOCUMENT gyotoy_export, filename

     Save current plot or Gyotoy data to file.

     The plot is exported using the standard Yorick functions if
     filename ends in ".pdf", ".png", ".eps", ".ps", ".jpeg", ".jpg",
     or ".jfif".

     The data are exported to an ASCII file if filename ends in ".dat"
     or ".txt". The GTK interface needs to be running for this to
     happen. Use gyotoy_save_data for the same result if the interface
     is not running.
     
   SEE ALSO: gyotoy_import, gyotoy_save_data
 */
  window, _gyotoy_wid;
  require,"pathfun.i";
  ext=pathsplit(filename,delim=".")(0);
  if (strglob("*.pdf", filename, case=0)) pdf, filename;
  else if (strglob("*.png", filename, case=0)) png, filename;
  else if (strglob("*.eps", filename, case=0)) eps, filename;
  else if (strglob("*.ps",  filename, case=0)) ps,  filename;
  else if (strglob("*.jpeg",filename, case=0)) jpeg,filename;
  else if (strglob("*.jpg", filename, case=0)) jpeg,filename;
  else if (strglob("*.jfif",filename, case=0)) jpeg,filename;
  else if (strglob("*.dat", filename, case=0)) gyotoy_save_data(filename);
  else if (strglob("*.txt", filename, case=0)) gyotoy_save_data(filename);
  else if (strglob("*.xml", filename, case=0)) _gyotoy_particle,
                                                 xmlwrite=filename;
  else pyk,"warning('Unknown file type for "+filename+"')";
}

func gyotoy_warning(msg) {
  pyk, "warning('"+msg+"')";
}

func gyotoy_import(filename) {
  extern _gyotoy_filename, _gyotoy_inhibit_redraw;
  // XML file:
  if (strpart(filename,-3:0)==".xml") {
    local rad2deg;
    rad2deg=180./pi;
    _gyotoy_filename=filename;
    gyotoy_set_particle, gyoto_Astrobj(filename);
    return;
  }

  // Non-xml file
  if(! (file = open(filename,"r",1))) {
    gyotoy_warning, "Could not open file "+filename+" for reading";
    return;
  }
  /*  if (!strmatch(rdline(file),"# Gyotoy save file")) {
    gyotoy_warning, filename+" is not a Gyotoy save file";
    return;
    }*/
  rdr=_gyotoy_inhibit_redraw;
  _gyotoy_inhibit_redraw=1;

  while ((line=rdline(file)) && !strmatch(line,"# Start Gyoto parameters")) {}
  if (!line)  {
    gyotoy_warning, filename+" is not a Gyoto save file";
    return;
  }
  while (!strmatch((line=rdline(file)),"# End Gyoto parameters")) {
    key="";
    value="";
    sread,line, format="# %s = %s", key, value;
    if (key=="metric_file") gyotoy_set_metric, strpart(value, 2:-1);
    ok=pyk("set_parameter('"+key+"',"+value+")");
  }
  _gyotoy_inhibit_redraw=rdr;
  gyotoy_compute_and_draw;
  //pyk,"compute_orbit('rien')";
}

func gyotoy_print(truc) {print,truc;}

func gyotoy_set_mass(mass) {
  extern _gyotoy_mass;
  _gyotoy_mass=mass;
  gyotoy_redraw;
}

func gyotoy_set_particle_type(type) {
  extern _gyotoy_particle_is_massless, _gyotoy_particle, _gyotoy_metric;
  extern _gyotoy_txyz;
  _gyotoy_particle_is_massless = (type=="photon");
  if (type=="photon") {
    _gyotoy_particle_is_massless = 1;
    _gyotoy_particle = gyoto_Photon();
  } else {
    _gyotoy_particle_is_massless = 0;
    _gyotoy_particle = gyoto_Star();
  }
  if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_particle,metric=_gyotoy_metric,
    initcoord=_gyotoy_initcoord(1:4),_gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_save_data(filename) {
/* DOCUMENT gyotoy_save_data,filename
   
     Save parameters and orbit to text file.
     
   SEE ALSO:
 */
  extern _gyotoy_metric;
  t=_gyotoy_txyz(,1);
  x=_gyotoy_txyz(,2);
  y=_gyotoy_txyz(,3);
  z=_gyotoy_txyz(,4);

  gyoto_convert, x, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
  gyoto_convert, y, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
  gyoto_convert, z, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
  ptype = (_gyotoy_particle_is_massless?"photon":"star");
  mtype = (_gyotoy_metric(kind=)=="KerrBL")?"kerrbl":"file";
  spin=(_gyotoy_metric(kind=)=="KerrBL")?_gyotoy_metric(spin=):2;
  if (f=open(filename,"w",1)) {
    write, f, format="# %s\n", "Gyoto save file";
    write, f, format="# %s\n", "Start Gyoto parameters";
    write, f, format="# %13s = \"%s\"\n", "particle_type", ptype;
    write, f, format="# %13s = \"%s\"\n", "metric_type", mtype;
    if (mtype=="kerrbl")
      write, f, format="# %13s = %14.12f\n", "spin", spin;
    else
      write, f, format="# %13s = \"%s\"\n", "metric_file", _gyotoy_metric_file;
    write, f, format="# %13s = %14.12f\n",
      ["mass",
       "t0", "r0", "theta0", "phi0",
       "rprime0", "thetaprime0", "phiprime0",
       "t1"],
      _([_gyotoy_mass],
        _gyotoy_initcoord,
        _gyotoy_t1);
    write, f, format="# %13s = \"%s\"\n", "length_unit", _gyotoy_unit;
    write, f, format="# %s\n", "End Gyoto parameters";
    write, f, format= "# %s\n", "Columns are t, x, y, z";
    write, f, format="%14.12f %14.12f %14.12f %14.12f\n", t, x, y, z;
    close,f;
  } else gyoto_warning, "Could not create file "+filename;
}

extern _pyk_proc, _gyotoy_stand_alone, _gyotoy_wid;
if (is_void(_gyotoy_wid)) _gyotoy_wid=0;
gyotoy_args=get_argv();
if (is_void(GYOTOY_NO_AUTO) & numberof(gyotoy_args)>=3 && anyof(basename(gyotoy_args(3))==["gyotoy.i","gyotoy"])) {
  _gyotoy_stand_alone=1;
  if (numberof(gyotoy_args)>3) {
    gyotoy_args=gyotoy_args(4:);
    ind=where(strgrep("^-",gyotoy_args)(2,)==-1);
    if (numberof(ind)) _gyotoy_filename=gyotoy_args(ind(1));
    if (numberof(ind)<numberof(gyotoy_args)) {
      ind=where(strgrep("^--",gyotoy_args)(2,)!=-1);
      gyotoy_options=gyotoy_args(ind);
      for (o=1;o<=numberof(gyotoy_options);o++) {
        option=gyotoy_options(o);
        pos=strfind("=",option);
        key=strpart(option,3:pos(1));
        if (pos(2)==-1) value= "true"; else value=strpart(option,pos(2)+1:);
        if (key == "stand-alone") {
          if (anyof(value==["true","1","TRUE","T","yes","t"])) _gyotoy_stand_alone=1;
          else _gyotoy_stand_alone=0;
        }
        else if (key == "ui") gyotoy_ui=value;
      }
    }
  }
  if (_gyotoy_stand_alone) batch,1;
  if (!is_void(_gyotoy_filename)) gyotoy,_gyotoy_filename;
  else gyotoy;
}
