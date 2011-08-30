#!/usr/local/bin/yorick -i
/*
    Copyright 2007 F. Rigaut
    Copyright 2011 Thibaut Paumard

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
  _gyotoy_t1=t1;
  if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_txyz=_gyotoy_particle(xfill=t1, get_txyz=1);
  gyotoy_redraw;
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
  _gyotoy_txyz=_gyotoy_particle(metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7),
                   xfill=_gyotoy_t1,
                   get_txyz=1);
  gyotoy_redraw;
}

func gyotoy_set_initcoord(t0, r0, theta0, phi0,
                          rprime0, thetaprime0, phiprime0) {
  extern _gyotoy_initcoord, _gyotoy_txyz;
  _gyotoy_initcoord=[t0, r0, theta0, phi0,
                     rprime0, thetaprime0, phiprime0];
  if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_txyz=_gyotoy_particle(initcoord=_gyotoy_initcoord(1:4),
                                _gyotoy_initcoord(5:7),
                                xfill=_gyotoy_t1,
                                get_txyz=1);
  gyotoy_redraw;
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

  extern _gyotoy_txyz;

  if (numberof(_gyotoy_txyz)>1) {
    x=_gyotoy_txyz(,2);
    y=_gyotoy_txyz(,3);
    z=_gyotoy_txyz(,4);
    gyoto_convert, x, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    gyoto_convert, y, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    gyoto_convert, z, _gyotoy_mass, _gyotoy_distance, _gyotoy_unit;
    rmax=sqrt(max(x^2+y^2+z^2));

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
  extern _gyotoy_wid, _gyotoy_wstyle, _gyotoy_filename, _gyotoy_particle;
  ok=pyk("sleep(0.1)")
  _gyotoy_wstyle="nobox.gs";
  window,_gyotoy_wid,wait=1,parent=parent_id,style="nobox.gs";
  limits, square=1;
  gnomon,1;
  cage3,1;
  //orient3,pi,-pi/2;
  pldefault, marks=0;
  if (!is_void(_gyotoy_particle_to_load)) gyotoy_set_particle,_gyotoy_particle_to_load;
  else if (_gyotoy_filename) gyotoy_import,_gyotoy_filename;
  else pyk,"compute_orbit('rien')";
  pyk,"orient3('rien')";
  pyk,"glade.get_widget('metric_type').set_active(0)";
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
  if (_gyotoy_stand_alone) quit;
}

func gyotoy(filename) {
/* DOCUMENT gyotoy [,filename]
         or gyotoy, star
         or gyotoy, photon
     Launch Gyotoy GTK interface
   SEE ALSO:
 */
  extern _pyk_proc, _gyotoy_filename, _gyotoy_running, _gyotoy_particle;
  extern _gyotoy_particle_to_load;
  _gyotoy_running=1;
  _gyotoy_filename=[];
  _gyotoy_particle_to_load=[];
  
  python_exec=find_in_path("gyotoy.py", takefirst=1,
                           path=pathform(_(get_cwd(),
                                           _(Y_SITES,
                                             Y_SITE)+"python/")));
 
  if (strpart(python_exec, 1:2)=="~/")
    python_exec=get_home()+"/"+strpart(python_exec, 3:);
  
  glade_file= find_in_path("gyotoy.glade", takefirst=1,
                           path=pathform(_("./", _(Y_SITES,Y_SITE)+"glade/")));

  if (strpart(glade_file, 1:2)=="~/")
    glade_file=get_home()+"/"+strpart(glade_file, 3:);
  
  gyotoytop = dirname(glade_file)+"/";
  
  // build command to spawn:
  pyk_cmd=[python_exec,gyotoytop];
  
  // spawn it and attach to _pyk_callback (see pyk.i):
  // that starts python and pass it the path to the glade file
  if (is_string(filename)) _gyotoy_filename=filename;
  else if (is_gyoto_Star(filename)||typeof(filename)=="gyoto_Photon")
    _gyotoy_particle_to_load=filename;
  else _gyotoy_particle_to_load=_gyotoy_particle;
  _pyk_proc = spawn(pyk_cmd, _pyk_callback);
  if (is_void(_pyk_proc)) error, "Failed to launch python";
  if (_gyotoy_filename &&
      (strglob("*.dat",_gyotoy_filename,case=0) |
       strglob("*.txt",_gyotoy_filename,case=0)))
    pyk,"set_filename('"+_gyotoy_filename+"')";
}

func gyotoy_set_particle(part) {
  extern _gyotoy_particle, _gyotoy_metric, _gyotoy_initcoord, _gyotoy_txyz;
  //  _gyotoy_particle=part;

  if (is_gyoto_Star(part)) {
    part_type="star";
    _gyotoy_particle_is_massless=0;
  }
  else {
    part_type="photon";
    _gyotoy_particle_is_massless=1;
  }

  ok=pyk("set_parameter('particle_type','"+part_type+"')");

  // Metric & projection
  _gyotoy_metric = part(metric=);
  if (_gyotoy_metric(kind=)=="KerrBL")
    ok=pyk("set_parameter('spin',"+swrite(format="%.12f",_gyotoy_metric(spin=))+")");
  //ok=pyk("set_parameter('incl',"+swrite(format="%.12f",metric(get_inclination=1)*rad2deg)+")");
  //ok=pyk("set_parameter('paln',"+swrite(format="%.12f",metric(get_paln=1)*rad2deg)+")");
  //ok=pyk("set_parameter('phase',"+swrite(format="%.12f",metric(get_argument=1)*rad2deg)+")");
  //ok=pyk("set_parameter('distance',"+swrite(format="%.12f",metric(get_distance=1))+")");
  //ok=pyk("set_parameter('mass',"+swrite(format="%.12f",metric(get_mass=1))+")");

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
  // XML file:
  if (strpart(filename,-3:0)==".xml") {
    local rad2deg;
    rad2deg=180./pi;
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
  while ((line=rdline(file)) && !strmatch(line,"# Start Gyoto parameters")) {}
  if (!line)  {
    gyotoy_warning, filename+" is not a Gyoto save file";
    return;
  }
  while (!strmatch((line=rdline(file)),"# End Gyoto parameters")) {
    key="";
    value="";
    sread,line, format="# %s = %s", key, value;
    ok=pyk("set_parameter('"+key+"',"+value+")");
  }
  pyk,"compute_orbit('rien')";
}

func gyotoy_print(truc) {print,truc;}

func gyotoy_set_mass(mass) {
  extern _gyotoy_mass;
  _gyotoy_mass=mass;
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
  _gyotoy_txyz=_gyotoy_particle(metric=_gyotoy_metric,
                                initcoord=_gyotoy_initcoord(1:4),
                                _gyotoy_initcoord(5:7),
                                xfill=_gyotoy_t1, get_txyz=1);
  gyotoy_redraw;
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
  spin=(_gyotoy_metric(kind=)=="KerrBL")?_gyotoy_metric(spin=):2;
  if (f=open(filename,"w",1)) {
    write, f, format="# %s\n", "Gyoto save file";
    write, f, format="# %s\n", "Start Gyoto parameters";
    write, f, format="# %13s = \"%s\"\n", "particle_type", ptype;
    write, f, format="# %13s = \"%s\"\n", "metric_type", "kerr";
    write, f, format="# %13s = %14.12f\n",
      ["spin", "mass",
       "t0", "r0", "theta0", "phi0",
       "rprime0", "thetaprime0", "phiprime0",
       "t1"],
      _([spin,
         _gyotoy_mass],
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
