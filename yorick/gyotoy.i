#!/usr/bin/env yorick -i
/*
    Copyright 2007 F. Rigaut
    Copyright 2011-2015 Thibaut Paumard

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

require, "gyoto.i";
require, "gyoto_std.i";
require, "pl3d.i";
// gy_gtk needs to be included after defining gy_gtk_before_init
// require, "gy_gtk.i";

func __gyotoy_before_init
{
  noop, gy.GLib.set_prgname("Gyotoy");
  noop, gy.GLib.set_application_name("Gyotoy");
  noop, gy.Gdk.set_allowed_backends("x11,*");
}

GYOTOY_VERSION = "0.2";

extern _gyotoy_running;
extern _gyotoy_reticle;
extern _gyotoy_wstyle;
extern _gyotoy_stand_alone;
extern _gyotoy_particle_is_massless;
extern _gyotoy_metric;
extern _gyotoy_particle;
extern _gyotoy_t1;
extern _gyotoy_delta;
extern _gyotoy_adaptive;
extern _gyotoy_initcoord;
extern _gyotoy_mass;
extern _gyotoy_distance;
extern _gyotoy_unit;
extern _gyotoy_txyz;
extern _gyotoy_particle_to_load;
extern _gyotoy_inhibit_redraw;
extern _gyotoy_metric_file;
extern _gyotoy_nsteps;
/* DOCUMENT extern _gyotoy_*;
     Some variable holding state information for use inside gyotoy.
     _running: true if GTK interface is running;
     _reticle: true if reticle should be displayed;
     .yid: Yorick window ID of the Gyotoy plot area;
     _wstyle: window style (normally "nobox.gs");
     _gyotoy.filename: name of the last file read (careful: can be
     the py file itself)
     _gyotoy_stand_alone: true if the user is not supposed to ever see
     the Yorick prompt.
     _gyotoy_particle_is_massless: true for a photon.
     _gyoto_metric: the actual metric object.
   SEE ALSO: gyotoy
 */
if (is_void(_gyotoy_nsteps)) _gyotoy_nsteps=100 ;

func gyotoy_adaptive(mode, void) {
  extern _gyotoy_adaptive;
  if (!is_numerical(mode)) mode = !mode.get_active();
  _gyotoy_adaptive=mode;
  if (mode) _gyotoy_particle, setparameter="Adaptive";
  else  _gyotoy_particle, setparameter="NonAdaptive";
}

func gyotoy_reset(void){
  extern _gyotoy_running;
  extern _gyotoy_reticle;
  extern _gyotoy_wstyle;
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
  _gyotoy_delta=0.01;
  _gyotoy_particle=gyoto_Star(metric=_gyotoy_metric,
                              initcoord=_gyotoy_initcoord(1:4),
                              _gyotoy_initcoord(5:7),
                              delta=_gyotoy_delta);
  gyotoy_adaptive,1;
  _gyotoy_particle, xfill=_gyotoy_t1;
  _gyotoy_txyz=_gyotoy_particle(get_txyz=1);
  _gyotoy_mass=4e6;
  _gyotoy_distance=8.;
  _gyotoy_unit="geometrical";
}
gyotoy_reset;

func gyotoy_set_unit(unit, data) {
  extern _gyotoy_unit;
  if (!is_string(unit)) {
    if (!unit.get_active()) return;
    unit=gy.Gtk.Buildable(unit).get_name();
  }
  _gyotoy_unit=unit;
  gyotoy_redraw;
}

func gyotoy_set_t1(t1, void)
{
  extern _gyotoy_txyz, _gyotoy_t1;
  if (!is_numerical(t1)) t1 = t1.get_value();
  /*  if ( (_gyotoy_t1 > _gyotoy_t0 && t1<_gyotoy_t1) ||
       (_gyotoy_t1 < _gyotoy_t0 && t1>_gyotoy_t1) )
       _gyotoy_particle, reset=;*/
  _gyotoy_t1=t1;
  if (catch(0x08)) return; // avoid breaking in case of v>c
  //  _gyotoy_txyz=_gyotoy_particle(xfill=t1, get_txyz=1);
  gyotoy_compute_and_draw;
}

func gyotoy_set_distance(dist, void) {
  extern _gyotoy_distance;
  if (!is_numerical(dist)) dist = dist.get_value();
  _gyotoy_distance=dist;
  gyotoy_redraw;
}

func gyotoy_set_KerrBL_metric(spin, void)
{
  extern _gyotoy_metric, _gyotoy_particle, _gyotoy_txyz;
  if (!is_numerical(spin)) spin = spin.get_value();

  _gyotoy_metric = gyoto_KerrBL ( spin = spin );
  if (_gyotoy_inhibit_redraw) return;
  //if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_particle,metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_set_KerrKS_metric(spin, void)
{
  extern _gyotoy_metric, _gyotoy_particle, _gyotoy_txyz;
  if (!is_numerical(spin)) spin = spin.get_value();

  _gyotoy_metric = gyoto_KerrKS ( spin = spin );
  if (_gyotoy_inhibit_redraw) return;
  //if (catch(0x08)) return; // avoid breaking in case of v>c
  _gyotoy_particle,metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_set_metric(fname, void)
{
  extern _gyotoy_metric, _gyotoy_particle, _gyotoy_txyz, _gyotoy_metric_file;

  if (catch(0x08)) {
    // avoid breaking in case bad file
    gyotoy_warning, "Unable to load metric. Is this a GYOTO XML description file?";
    return;
  }
  if (is_void(structof(fname)))
    fname = gy.Gtk.FileChooser(fname).get_filename();
  if (is_string(fname)) metric=gyoto_Metric(fname);
  if (typeof(metric)!="gyoto_Metric") gyotoy_warning, "Failed to set metric";
  prev_metric=_gyotoy_metric;
  _gyotoy_metric=metric;
  _gyotoy_metric_file=fname;
  if (catch(0x08)) {
    // avoid breaking in case of v>c or other problem
    gyotoy_warning, "metric loaded but orbit computation failed";
    _gyotoy_metric=prev_metric;
    _gyotoy_particle,metric=_gyotoy_metric,
      initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
    return;
  }
  _gyotoy_particle,metric=_gyotoy_metric,
                   initcoord=_gyotoy_initcoord(1:4), _gyotoy_initcoord(5:7);
  gyotoy_compute_and_draw;
}

func gyotoy_set_initcoord(t0, r0, theta0, phi0,
                          rprime0, thetaprime0, phiprime0) {
  extern _gyotoy_initcoord, _gyotoy_txyz, _gyotoy_inhibit_redraw;
  if (!is_numerical(t0)) {
    r0          = _gyotoy.builder.get_object("r0").get_value();
    theta0      = _gyotoy.builder.get_object("theta0").get_value();
    phi0        = _gyotoy.builder.get_object("phi0").get_value();
    t0          = _gyotoy.builder.get_object("t0").get_value();
    rprime0     = _gyotoy.builder.get_object("rprime0").get_value();
    thetaprime0 = _gyotoy.builder.get_object("thetaprime0").get_value();
    phiprime0   = _gyotoy.builder.get_object("phiprime0").get_value();
  }
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

    window, _gyotoy.yid;
    clear3;
    if (_gyotoy_reticle) gyoto_reticle, rmax*1.41421;
    if (gyotoy_origin) gyoto_plmksky,0,0;
    gyoto_plg3, x, y, z;
    limit3,-rmax,rmax,-rmax,rmax,-rmax,rmax;
    //if (is_void(set_idler())) draw3,1;
    //gy_gtk_idleonce;
  }

}

func gyotoy_toggle_reticle(wdg, data) {
/* DOCUMENT gyotoy_toggle_reticle
    Toggles sky reticle on/off
   SEE ALSO: gyoto_reticle
 */
  extern _gyotoy_reticle;
  if (is_void(_gyotoy_reticle)) _gyotoy_reticle=1;
  else _gyotoy_reticle=1-_gyotoy_reticle;
  gyotoy_redraw;
}

func gyotoy_window_init
// initialize window in GTK frontend
{
  extern _gyotoy_wstyle, _gyotoy_particle, _gyotoy_inhibit_redraw;

  iconf = find_in_path("gyotoy.png", takefirst=1,
                       path=pathform(_(get_cwd(),
                                       _(Y_SITES,
                                         Y_SITE)+"data/")));
  
  if (iconf) {
    icon = gy.GdkPixbuf.Pixbuf.new_from_file(iconf);
    noop, _gyotoy.toplevel.set_icon(icon);
  }

  limits, square=1;
  gnomon,1;
  cage3,1;
  //orient3,pi,-pi/2;
  pldefault, marks=0;
  if (!is_void(_gyotoy_particle_to_load)) gyotoy_set_particle,_gyotoy_particle_to_load;
  else if (_gyotoy.filename) gyotoy_import,_gyotoy.filename;
  //else pyk,"compute_orbit('rien')";
  _gyotoy_orient3;
  noop, _gyotoy.builder.get_object("metric_type").set_active(0);
  _gyotoy_inhibit_redraw=0;
  gyotoy_compute_and_draw;
}

func gyotoy_quit(wdg, void1, void2) {
// called when GTK window is closed
  extern _gyotoy_running;
  if (_gyotoy_stand_alone) quit;
  noop, _gyotoy.toplevel.hide();
  _gyotoy_running=0;
  gy_gtk_idler_maybe_stop;
}

func gyotoy_checkvers_cb(msg) {
  extern _gyotoy_pyvers;
  _gyotoy_pyvers=msg;
}

if (is_void(_gyotoy)) _gyotoy = save();

func gyotoy(filename) {
/* DOCUMENT gyotoy [,filename]
         or gyotoy, star
         or gyotoy, photon
     Launch Gyotoy GTK interface
   SEE ALSO:
 */
  require, "gy_gtk.i";
  extern _gyotoy, _gyotoy_running, _gyotoy_particle;
  extern _gyotoy_particle_to_load, _gyotoy_inhibit_redraw;
  local yid;
  yid=[];

  _gyotoy=save();
  _gyotoy_running=1;
  _gyotoy, filename=[];
  _gyotoy_particle_to_load=[];
  _gyotoy_inhibit_redraw=1;

  _gyotoy, builder = gy_gtk_builder("gyotoy.xml");

  if (is_string(filename)) _gyotoy, filename=filename;
  else if ((is_gyoto_Astrobj(filename) && filename(kind=)=="Star")
           ||typeof(filename)=="gyoto_Photon")
    _gyotoy_particle_to_load=filename;
  else _gyotoy_particle_to_load=_gyotoy_particle;

  _gyotoy, toplevel = _gyotoy.builder.get_object("main_window");
  noop, _gyotoy.builder.get_object("plot_slot").
    add(gy_gtk_ywindow(yid, style="nobox.gs",on_realize=gyotoy_window_init));
  _gyotoy, yid = yid;
  noop, _gyotoy.builder.get_object("main_vbox").
    pack_start(gy_gtk_ycmd(), 0, 0, 0);
  gy_signal_connect, _gyotoy.builder;
  gy_signal_connect, _gyotoy.toplevel, "delete-event", gyotoy_quit;

  // git has gy_gtk_main, but it's not yet in the stable release
  if (is_func(gy_gtk_main)) gy_gtk_main, _gyotoy.toplevel;
  else noop, _gyotoy.toplevel.show_all();
  
  /*
  if (_gyotoy.filename &&
      (strglob("*.dat",_gyotoy.filename,case=0) ||
       strglob("*.txt",_gyotoy.filename,case=0) ||
       strglob("*.xml",_gyotoy.filename,case=0) ))
    pyk,"set_filename('"+_gyotoy.filename+"')";
  */
}

func _gyotoy_on_realize
{
  extern _gyotoy_inhibit_redraw;
  _gyotoy_inhibit_redraw=0;
  gyotoy_redraw;
}

func gyotoy_set_particle(part) {
  tp = typeof(part);
  if (!((tp == "gyoto_Astrobj" && part(kind=) == "Star") ||
        (tp == "gyoto_Photon"))) {
    error, "Particle must be Star or Photon";
    return;
  }
  
  extern _gyotoy_particle, _gyotoy_metric, _gyotoy_initcoord, _gyotoy_txyz,
    _gyotoy_metric_file, _gyotoy_inhibit_redraw,
    _gyotoy_delta, _gyotoy_adaptive;
  rdr=_gyotoy_inhibit_redraw;
  _gyotoy_inhibit_redraw=1;
  
  omtype=_gyotoy_metric(kind=);
  oldmass=_gyotoy_particle_is_massless;
  
  _gyotoy_particle=part;
  _gyotoy_metric = part(metric=);
  _gyotoy_initcoord=part(initcoord=);
  _gyotoy_delta=part(delta=);
  _gyotoy_adaptive=part(adaptive=);
  
  if (is_gyoto_Astrobj(part)) {
    part_type="star";
    _gyotoy_particle_is_massless=0;
  }
  else {
    part_type="photon";
    _gyotoy_particle_is_massless=1;
  }
  if (oldmass != _gyotoy_particle_is_massless)
    _noop, gyotoy.builder.get_object(part_type).set_active(1);

  noop, _gyotoy.builder.get_object("delta").set_value(_gyotoy_delta);
  noop, _gyotoy.builder.get_object("adaptive_button").
    set_active(1-_gyotoy_adaptive);

  // Metric & projection
  if ((mtype=_gyotoy_metric(kind=))=="KerrBL") {
    if (omtype != "KerrBL")
      noop, _gyotoy.builder.get_object("metric_type").set_active(0);
    noop, _gyotoy.builder.get_object("spin").set_value(_gyotoy_metric(spin=));
  } else {
    if (otype == "KerrBL")
      noop, _gyotoy.builder.get_object("metric_type").set_active(1);
    if (_gyotoy.filename)
      noop, gy.Gtk.FileChooser(_gyotoy.builder.get_object("metric_file")).
        set_filename(_gyotoy.filename);
    _gyotoy_metric_file=_gyotoy.filename;
  }
  
  //ok=pyk("set_parameter('incl',"+swrite(format="%.12f",metric(get_inclination=1)*rad2deg)+")");
  //ok=pyk("set_parameter('paln',"+swrite(format="%.12f",metric(get_paln=1)*rad2deg)+")");
  //ok=pyk("set_parameter('phase',"+swrite(format="%.12f",metric(get_argument=1)*rad2deg)+")");
  //ok=pyk("set_parameter('distance',"+swrite(format="%.12f",metric(get_distance=1))+")");
  //  m_sun = 1.98843e30;     // kg
  //  ok=pyk("set_parameter('mass',"+swrite(format="%.12f",_gyotoy_metric(mass=)/m_sun)+")");

  // Initial condition
  coord = part(initcoord=);
  noop, _gyotoy.builder.get_object("t0").set_value(coord(1));
  noop, _gyotoy.builder.get_object("r0").set_value(coord(2));
  noop, _gyotoy.builder.get_object("theta0").set_value(coord(3));
  noop, _gyotoy.builder.get_object("phi0").set_value(coord(4));
  
  if (part_type=="star") fact=1./coord(5);
  else fact=1.;
  noop, _gyotoy.builder.get_object("rprime0").set_value(coord(6)*fact);
  noop, _gyotoy.builder.get_object("thetaprime0").set_value(coord(7)*fact);
  noop, _gyotoy.builder.get_object("phiprime0").set_value(coord(8)*fact);

  // Wait for parameters to have reached glade
  //    pyk,"compute_orbit('rien')";

  // bug ?
  if (mtype != "KerrBL")
    noop, _gyotoy.builder.get_object("metric_type").set_active(1);


  _gyotoy_inhibit_redraw=rdr;
  gyotoy_compute_and_draw;
  
}

func gyotoy_set_nsteps(nsteps, void) {
  extern _gyotoy_nsteps;
  if (is_numerical(nsteps)) nsteps=long(nsteps);
  else nsteps = nsteps.get_value();
  if (nsteps <=0) nsteps=1;
  _gyotoy_nsteps=nsteps;
}

func gyotoy_compute_and_draw(rien) {
  
  extern _gyotoy_particle, _gyotoy_redrawing, _gyotoy_cancel, _gyotoy_nsteps;
  extern _gyotoy_t1, _gyotoy_inhibit_redraw, _gyotoy_txyz;
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
    gyotoy_set_play_image, "gtk-media-play";
    gyotoy_set_fraction, (t-t0)/(_gyotoy_t1-t0);
    gyotoy_redraw;
    --_gyotoy_redrawing;
    return;
  }

  gyotoy_set_fraction, (t-t0)/(_gyotoy_t1-t0);
  gyotoy_set_play_image, "gtk-media-pause";
  
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
    //gy_gtk_step_once;
    frac=(t-t0)/(_gyotoy_t1-t0);
    if (frac>1.) frac=1;
    gyotoy_set_fraction, frac;
    while (gy.Gtk.events_pending ())  noop, gy.Gtk.main_iteration ();
  }

  gyotoy_set_play_image, "gtk-media-play";
  frac=(t-t0)/(_gyotoy_t1-t0);
  if (frac>1.) frac=1;
  gyotoy_set_fraction, frac;

  gyotoy_redraw;
  --_gyotoy_redrawing;
  if (_gyotoy_cancel==2) gyotoy_rewind;
}

func gyotoy_set_play_image(name)
{
  noop, _gyotoy.builder.get_object("play_image").
    set_from_stock(name, gy.Gtk.IconSize.button);
}

func gyotoy_rewind(wdg, data)
{
  extern _gyotoy_cancel;
  if (_gyotoy_redrawing) {
    _gyotoy_cancel=2;
    return;
  }
  _gyotoy_particle, reset=;
  gyotoy_set_fraction, 0;
}

func gyotoy_set_fraction(frac)
{
  noop, _gyotoy.builder.get_object("progressbar").set_fraction(frac);
}                          

func gyotoy_set_delta(delta, void) {
  extern _gyotoy_delta;
  if (!is_numerical(delta)) delta = delta.get_value();
  _gyotoy_delta=delta;
  _gyotoy_particle, delta=delta;
}

func gyotoy_inhibit_redraw(mode, data) {
  extern _gyotoy_inhibit_redraw, _gyotoy_cancel;
  if (!structof(mode)) mode=mode.get_active();
  _gyotoy_inhibit_redraw=mode;
  if (mode &&_gyotoy_redrawing) _gyotoy_cancel=1;
}

func gyotoy_play_pause(wdg, data)
{
  noop, _gyotoy.builder.get_object("inhibit_button").set_active(0);
  gyotoy_compute_and_draw;
}

func gyotoy_export(filename, data)
{
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

  // If filename is not a string, assume open button has been clicked
  if (!is_string(filename)) {
    Gtk = gy.require("Gtk", "3.0");

    chooser = Gtk.FileChooserDialog();
    noop, chooser.add_button(Gtk.STOCK_CANCEL, Gtk.ResponseType.cancel);
    noop, chooser.add_button(Gtk.STOCK_SAVE, Gtk.ResponseType.ok);
    fcfc = Gtk.FileChooser(chooser);
    noop, fcfc.set_action(Gtk.FileChooserAction.save);
    noop, fcfc.set_do_overwrite_confirmation(1);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*");
    noop, filter.set_name("All Files");
    noop, fcfc.add_filter(filter);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*.[tT][xX][tT]");
    noop, filter.add_pattern("*.[dD][aA][tT]");
    noop, filter.add_pattern("*.[pP][dD][fF]");
    noop, filter.add_pattern("*.[eE][pP][sS]");
    noop, filter.add_pattern("*.[jJ][pP][eE][gG]");
    noop, filter.add_pattern("*.[jJ][pP][gG]");
    noop, filter.add_pattern("*.[jJ][fF][iI][fF]");
    noop, filter.add_pattern("*.[pP][nN][gG]");
    noop, filter.set_name("All supported files");
    noop, fcfc.add_filter(filter);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*.[dD][aA][tT]");
    noop, filter.add_pattern("*.[tT][xX][tT]");
    noop, filter.set_name("Data Files");
    noop, fcfc.add_filter(filter);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*.[xX][mM][lL]");
    noop, filter.set_name("XML description files");
    noop, fcfc.add_filter(filter);

   filter = Gtk.FileFilter();
   noop, filter.add_pattern("*.[pP][dD][fF]");
   noop, filter.set_name("PDF documents");
   noop, fcfc.add_filter(filter);

   filter = Gtk.FileFilter();
   noop, filter.add_pattern("*.[eE][pP][sS]");
   noop, filter.set_name("EPS documents");
   noop, fcfc.add_filter(filter);

   filter = Gtk.FileFilter();
   noop, filter.add_pattern("*.[jJ][pP][eE][gG]");
   noop, filter.add_pattern("*.[jJ][pP][gG]");
   noop, filter.add_pattern("*.[jJ][fF][iI][fF]");
   noop, filter.set_name("JPEG images");
   noop, fcfc.add_filter(filter);

   filter = Gtk.FileFilter();
   noop, filter.add_pattern("*.[pP][nN][gG]");
   noop, filter.set_name("PNG images");
   noop, fcfc.add_filter(filter);

   res = chooser.run();
    noop, chooser.hide();
    if (res == Gtk.ResponseType.OK) {
      filename=fcfc.get_filename();
      noop, chooser.destroy();
    } else {
      noop, chooser.destroy();
      return;
    }
  }

  window, _gyotoy.yid;
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

func gyotoy_import(filename, data) {
  extern _gyotoy, _gyotoy_inhibit_redraw;

  // If filename is not a string, assume open button has been clicked
  if (!is_string(filename)) {
    Gtk = gy.require("Gtk", "3.0");

    chooser = Gtk.FileChooserDialog();
    noop, chooser.add_button(Gtk.STOCK_CANCEL, Gtk.ResponseType.cancel);
    noop, chooser.add_button(Gtk.STOCK_OPEN, Gtk.ResponseType.ok);
    fcfc = Gtk.FileChooser(chooser);
    noop, fcfc.set_action(Gtk.FileChooserAction.open);
    noop, fcfc.set_do_overwrite_confirmation(1);
    /*
    buttons=(Gtk.STOCK_CANCEL,
                                     Gtk.ResponseType.CANCEL,
                                     Gtk.STOCK_OPEN,
                                     Gtk.ResponseType.OK));
    */
    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*");
    noop, filter.set_name("All Files");
    noop, fcfc.add_filter(filter);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*.[dD][aA][tT]");
    noop, filter.add_pattern("*.[tT][xX][tT]");
    noop, filter.set_name("Data Files");
    noop, fcfc.add_filter(filter);

    filter = Gtk.FileFilter();
    noop, filter.add_pattern("*.[xX][mM][lL]");
    noop, filter.set_name("XML description files");
    noop, fcfc.add_filter(filter);

    res = chooser.run();
    noop, chooser.hide();
    if (res == Gtk.ResponseType.OK) {
      filename=fcfc.get_filename();
      noop, chooser.destroy();
    } else {
      noop, chooser.destroy();
      return;
    }
  }
           
  // XML file:
  if (strpart(filename,-3:0)==".xml") {
    local rad2deg;
    rad2deg=180./pi;
    _gyotoy, filename=filename;
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

func gyotoy_set_mass(mass, void)
{
  extern _gyotoy_mass;
  if (!is_numerical(mass)) mass = mass.get_value();
  _gyotoy_mass=mass;
  gyotoy_redraw;
}

func gyotoy_set_particle_type(type, void)
{
  extern _gyotoy_particle_is_massless, _gyotoy_particle, _gyotoy_metric;
  extern _gyotoy_txyz;
  if (!is_string(type)) {
    if (!type.get_active()) return;
    type = gy.Gtk.Buildable(type).get_name();
  }
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
  _gyotoy_particle, delta=_gyotoy_delta;
  gyotoy_adaptive, _gyotoy_adaptive;
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
       "t1", "delta"],
      _([_gyotoy_mass],
        _gyotoy_initcoord,
        _gyotoy_t1, _gyotoy_delta);
    write, f, format="# %13s = %d\n",
      ["adaptive", "nsteps"], [_gyotoy_adaptive, _gyotoy_nsteps];
    write, f, format="# %13s = \"%s\"\n", "length_unit", _gyotoy_unit;
    write, f, format="# %s\n", "End Gyoto parameters";
    write, f, format= "# %s\n", "Columns are t, x, y, z";
    write, f, format="%14.12f %14.12f %14.12f %14.12f\n", t, x, y, z;
    close,f;
  } else error, "Could not create file "+filename;
}


/************** CALLBACKS ****************/
func gyotoy_about(wdg, udata)
{
  Gtk=gy.require("Gtk", "3.0");
  dialog = Gtk.AboutDialog();
  noop, dialog.set_program_name("Gyotoy");
  noop, dialog.set_version(GYOTOY_VERSION);
  noop, dialog.set_logo(icon);
  noop, dialog.set_copyright("Copyright © 20011-2013 Frédéric Vincent & Thibaut Paumard");
  noop, dialog.set_license_type(Gtk.License.gpl_3_0);
  noop, dialog.run();
  noop, dialog.destroy();
}

func _gyotoy_orient3(wdg, data)
{
  incl = _gyotoy.builder.get_object("incl").get_value();
  paln = _gyotoy.builder.get_object("paln").get_value();
  phase = _gyotoy.builder.get_object("phase").get_value();
  gyoto_orient3, incl, paln, phase;
  //draw3,1;
}

func _gyotoy_metric_type(wdg, data)
{
  type_id = wdg.get_active();
  if (type_id == 0) type = "kerrbl";
  else if (type_id == 1) type = "kerrks";
  else if (type_id == 2) type = "file";

  if (type=="file") {
    noop, _gyotoy.builder.get_object("spin").set_sensitive(0);
    noop, _gyotoy.builder.get_object("spin").set_visible(0);
    noop, _gyotoy.builder.get_object("spin_label").set_visible(0);
    noop, _gyotoy.builder.get_object("metric_file").set_sensitive(1);
    noop, _gyotoy.builder.get_object("metric_file_label").set_visible(1);
    noop, _gyotoy.builder.get_object("metric_file").set_visible(1);
    // TO PORT
    //_gyotoy.metric_file_set_cb(_gyotoy.builder.get_object("metric_file"));
  } else if (type=="kerrbl") {
    noop, _gyotoy.builder.get_object("spin").set_visible(1);
    noop, _gyotoy.builder.get_object("spin").set_sensitive(1);
    noop, _gyotoy.builder.get_object("spin_label").set_visible(1);
    noop, _gyotoy.builder.get_object("metric_file").set_sensitive(0);
    noop, _gyotoy.builder.get_object("metric_file_label").set_visible(0);
    noop, _gyotoy.builder.get_object("metric_file").set_visible(0);
    gyotoy_set_KerrBL_metric, _gyotoy.builder.get_object("spin");
  } else if (type=="kerrks") {
    noop, _gyotoy.builder.get_object("spin").set_visible(1);
    noop, _gyotoy.builder.get_object("spin").set_sensitive(1);
    noop, _gyotoy.builder.get_object("spin_label").set_visible(1);
    noop, _gyotoy.builder.get_object("metric_file").set_sensitive(0);
    noop, _gyotoy.builder.get_object("metric_file_label").set_visible(0);
    noop, _gyotoy.builder.get_object("metric_file").set_visible(0);
    gyotoy_set_KerrKS_metric, _gyotoy.builder.get_object("spin");
  } else {
    error, "Unknown metric type";
  }
}

func _gyotoy_gnomon(wdg, data)
{
  gnomon;
}

func _gyotoy_cage(wdg, data)
{
  cage3;
}

func _gyotoy_unlimit(wdg, data)
{
  limits;
}


extern gyotoy_warning;
gyotoy_warning=error;

extern _pyk_proc, _gyotoy_stand_alone;
if (is_void(_gyotoy.yid)) _gyotoy,yid=0;
gyotoy_args=get_argv();
if (is_void(GYOTOY_NO_AUTO) & numberof(gyotoy_args)>=3 && anyof(basename(gyotoy_args(3))==["gyotoy.i","gyotoy"])) {
  _gyotoy_stand_alone=1;
  if (numberof(gyotoy_args)>3) {
    gyotoy_args=gyotoy_args(4:);
    ind=where(strgrep("^-",gyotoy_args)(2,)==-1);
    if (numberof(ind)) _gyotoy, filename=gyotoy_args(ind(1));
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
  
  if (_gyotoy_stand_alone) {
    batch,1;
    gy_gtk_before_init = __gyotoy_before_init;
  }
  
  if (!is_void(_gyotoy.filename)) gyotoy,_gyotoy.filename;
  else gyotoy;
}
