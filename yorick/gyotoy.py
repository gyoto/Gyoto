#!/usr/bin/env python
# authors: Thibaut Paumard
# gyotoy.py
# Default parameters:
# Start Gyoto parameters
# particle_type = "star"
#   metric_type = 0
#          spin = 0.995000000000
#            t0 = 0.000000000000
#            r0 = 10.791000000000
#        theta0 = 1.570800000000
#          phi0 = 0.000000000000
#       rprime0 = 0.000000000000
#   thetaprime0 = 0.000000000000
#     phiprime0 = 0.016664000000
#            t1 = 3000.000000000000
#          incl = 120.000000000000
#          paln = -180.000000000000
#         phase = -120.000000000000
# End Gyoto parameters


#  Copyright 2007 F. Rigaut
#  Copyright 2011 Thibaut Paumard
#
#  This file is part of Gyoto.
#
#  Gyoto is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Gyoto is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.



# import necessary modules
import gtk
import gtk.glade
import sys
import gobject
import os, fcntl, errno
import time

# define the class
class gyotoy:
   
   def destroy(self, wdg, data=None):
      self.py2yo('gyotoy_quit')
      gtk.main_quit() # we'll call this from yorick before it quits.
      
   def __init__(self,gyotoytop):
      gobject.set_application_name("Gyotoy")
      gobject.set_prgname("Gyotoy")

      self.gyotoytop = gyotoytop
      self._pyk_blocked=0
      self.length_unit="geometrical"

      # read GUI definition
      self.glade = gtk.glade.XML(os.path.join(self.gyotoytop,'gyotoy.glade')) 
      
      # handle destroy event (so that we can kill the window through window bar)
      self.window = self.glade.get_widget('window1')
      if (self.window):
         self.window.connect('destroy', self.destroy)
         
      # autoconnect to callback functions
      # this will automatically connect the event handler "on_something_event"
      # to the here-defined function of the same name. This avoid defining a
      # long dictionary that serves the same purpose
      self.glade.signal_autoconnect(self)

      # set stdin non blocking, this will prevent readline to block
      # stdin is coming from yorick (yorick spawned this python process)
      fd = sys.stdin.fileno()
      flags = fcntl.fcntl(fd, fcntl.F_GETFL)
      fcntl.fcntl(fd, fcntl.F_SETFL, flags | os.O_NONBLOCK)
      
      # ... and add stdin to the event loop (yorick input pipe by spawn)
      gobject.io_add_watch(sys.stdin,gobject.IO_IN|gobject.IO_HUP,self.yo2py,None)

      # run: realize the interface, start event management
      gtk.main()

   # callback functions
   # This is the important, user-defined heart of this file.
   # here, you define the function that will be called when an event is raised
   # (read: happens): a button is pressed, a slider is moved, you enter a
   # drawingarea or move the cursor, etc... Just define a handler (callback) of
   # the same name in the glade UI definition (see signals) and you're done.

   def set_metric(self, wdg):
      spin        = self.glade.get_widget('spin').get_value()
      self.py2yo('gyotoy_set_KerrBL_metric %14.12f' % (spin))

   def set_mass(self, wdg):
      mass        = self.glade.get_widget('mass').get_value()
      self.py2yo('gyotoy_set_mass %14.12f' % (mass))

   def set_initcoord(self, wdg):
      r0          = self.glade.get_widget('r0').get_value()
      theta0      = self.glade.get_widget('theta0').get_value()
      phi0        = self.glade.get_widget('phi0').get_value()
      t0          = self.glade.get_widget('t0').get_value()
      rprime0     = self.glade.get_widget('rprime0').get_value()
      thetaprime0 = self.glade.get_widget('thetaprime0').get_value()
      phiprime0   = self.glade.get_widget('phiprime0').get_value()
      self.py2yo('gyotoy_set_initcoord %14.12f %14.12f %14.12f %14.12f %14.12f %14.12f %14.12f' % (t0, r0, theta0, phi0, rprime0, thetaprime0, phiprime0))

   def set_distance(self,wdg):
      distance    = self.glade.get_widget('distance').get_value()
      self.py2yo('gyotoy_set_distance %14.12f' % (distance))
      
   def set_t1(self,wdg):
      t1          = self.glade.get_widget('t1').get_value()
      self.py2yo('gyotoy_set_t1 %14.12f' % (t1))
      
   def orient3(self, wdg):
      incl = self.glade.get_widget('incl').get_value()
      paln = self.glade.get_widget('paln').get_value()
      phase = self.glade.get_widget('phase').get_value()
      self.py2yo("gyoto_orient3 %14.12f %14.12f %14.12f" % (incl, paln, phase))

   def sleep(self,delay):
      time.sleep(delay)
      self.pyk_resume("1")

   def metric_type_map(self,wdg, *args):
      self.glade.get_widget('metric_type').set_active(0)

   def on_drawingarea_map_event(self,wdg,*args):
#      time.sleep(2)
      self.glade.get_widget('metric_type').set_active(0)
      drawingarea = self.glade.get_widget('yorick_window')
      mwid = drawingarea.window.xid;
      self.py2yo('gyotoy_window_init %d' % mwid)
      self.orient3(wdg)
      #self.outfile=""
      # optionaly (but this is a good place to
      # do it), update GUI parameters from yorick:
      # self.py2yo('gui_update')     
      # minimal wrapper for yorick/python communication
      # This is really internal to the yorick/python communication and you should
      # not have to fiddle with it

   def on_save_as_activate(self,wdg):
      chooser = gtk.FileChooserDialog(title="Save plot or data to file",action=gtk.FILE_CHOOSER_ACTION_SAVE,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))

      filter = gtk.FileFilter()
      filter.add_pattern('*')
      filter.set_name('All Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[pP][dD][fF]')
      filter.add_pattern('*.[pP][nN][gG]')
      filter.add_pattern('*.[jJ][pP][gG]')
      filter.add_pattern('*.[eE][pP][sS]')
      filter.set_name('Supported image files (PDF, PNG, JPG and EPS)')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[dD][aA][tT]')
      filter.add_pattern('*.[tT][xX][tT]')
      filter.set_name('Data Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[pP][dD][fF]')
      filter.set_name('PDF Image Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[pP][nN][gG]')
      filter.set_name('PNG Image Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[jJ][pP][gG]')
      filter.set_name('JPG Image Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[eE][pP][sS]')
      filter.set_name('EPS Image Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[xX][mM][mL]')
      filter.set_name('XML description files')
      chooser.add_filter(filter)

      res = chooser.run()
      if res == gtk.RESPONSE_OK:
         file=chooser.get_filename()
         self.set_filename(file)
         self.py2yo('gyotoy_export "'+file+'"')
      chooser.destroy()

   def on_save_activate(self,wdg):
      self.py2yo('gyotoy_export "'+self.outfile+'"')

   def on_open_activate(self,wdg):
      chooser = gtk.FileChooserDialog(title="Open Gyotoy data file",action=gtk.FILE_CHOOSER_ACTION_OPEN,buttons=(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL,gtk.STOCK_OPEN,gtk.RESPONSE_OK))

      filter = gtk.FileFilter()
      filter.add_pattern('*')
      filter.set_name('All Files')
      chooser.add_filter(filter)

      filter = gtk.FileFilter()
      filter.add_pattern('*.[dD][aA][tT]')
      filter.add_pattern('*.[tT][xX][tT]')
      filter.set_name('Data Files')
      chooser.add_filter(filter)

      res = chooser.run()
      if res == gtk.RESPONSE_OK:
         file=chooser.get_filename()
         self.set_filename(file)
         self.py2yo('gyotoy_import "'+file+'"')
      chooser.destroy()

   def set_parameters_text(self, parameters, values):
      self.set_parameter_text(parameters[0],values[0])
      self.set_parameter_text(parameters[1],values[1])
      self.set_parameter_text(parameters[2],values[2])
      self.set_parameter_text(parameters[3],values[3])
      self.set_parameter_text(parameters[4],values[4])

   def set_length_unit(self, wdg):
      if wdg.get_active():
         self.py2yo('gyotoy_set_unit "%s"' % (wdg.get_name()))

   def set_particle_type(self, wdg):
      self.py2yo('gyotoy_set_particle_type "%s" ' % (wdg.get_name()))

   def set_parameter_text(self, param, value):
      self.glade.get_widget(param).set_text(value)

   def set_parameter(self, param, value):
      if (param=='metric_type'):
         self.glade.get_widget('metric_type').set_active(0)         
      elif (param=='length_unit' or param=='particle_type'):
         self.glade.get_widget(value).set_active(1)
      else:
         self.glade.get_widget(param).set_value(value)
         count=0
         while (self.glade.get_widget(param).get_value()!=value):
            time.sleep(1)
            count+=1
            if (count==5):
               self.py2yo('print "%s" %.12f %.12f' % (param, value, self.glade.get_widget(param).get_value()))
               self.yerror('Unable to set value')
      if (param=='spin'):
         self.set_metric(self.glade.get_widget(param))
      self.pyk_resume('1')

   def toggle_gnomon(self,wdg):
      self.py2yo('gnomon')

   def toggle_cage(self,wdg):
      self.py2yo('cage3')

   def limits(self,wdg):
      self.py2yo('limits')

   def toggle_reticle(self,wdg):
      self.py2yo('gyotoy_toggle_reticle')
      #mwid = self.glade.get_widget('yorick_window').window.xid;
      #self.py2yo('gyotoy_toggle_window_style %d' % mwid)
      self.redraw('rien')

   def redraw(self,wdg):
      self.py2yo('gyotoy_redraw')


   def about_window(self,wdg):
      dialog = self.glade.get_widget('aboutdialog')
      dialog.run()
      dialog.hide()

   def warning(self,msg):
      mbox = gtk.MessageDialog(self.window, gtk.DIALOG_MODAL, gtk.MESSAGE_WARNING, gtk.BUTTONS_OK, 'Gyotoy Warning');
      mbox.format_secondary_markup(msg);
      #,message-type=gtk.MESSAGE_WARNING,buttons=gtk.BUTTONS_OK);
      res=mbox.run();
      mbox.destroy();

   def set_filename(self,file):
      self.outfile=file
      self.glade.get_widget('save_button').set_sensitive(1)
      self.window.set_title("Gyotoy - "+file)

# pyk functions

   def pyk_sync(self):
      self._pyk_blocked=1
      sys.stdout.write('-s+y-n+c-+p-y+k-\n')
      sys.stdout.flush()

   def pyk_resume(self,msg):
      if (self._pyk_blocked):
         self.py2yo('pyk_resume ' + msg)
         self._pyk_blocked=0

   def py2yo(self,msg):
      # sends string command to yorick
      sys.stdout.write(msg+'\n')
      sys.stdout.flush()
   
   def yo2py(self,cb_condition,*args):
      if cb_condition == gobject.IO_HUP:
         raise SystemExit, "lost pipe to yorick"
      # handles string command from yorick
      # note: inidividual message needs to end with \n for proper ungarbling
      while 1:
         try:
            msg = sys.stdin.readline()
            msg = "self."+msg
            exec(msg)
         except IOError, e:
            if e.errno == errno.EAGAIN:
               # the pipe's empty, good
               break
            # else bomb out
            self.yerror(str(e))
            raise SystemExit, "yo2py unexpected IOError:" + str(e)
         except Exception, e:
#            raise SystemExit, "yo2py unexpected Exception:" + str(ee)
            self.yerror(str(e))
            raise SystemExit, "yo2py unexpected Exception:" + str(ee) +msg
         return True

   def yerror(self, msg):
         if (self._pyk_blocked):
            self.pyk_resume('0')
         self.py2yo('error "%s"' % msg)
         
# check calling syntax (should not be a problem as it is always called
# from yorick
      
if len(sys.argv) != 2:
   print 'Usage: gyotoy.py path_to_glade'
   raise SystemExit

# get path to glade file
gyotoytop = str(sys.argv[1])
# execute it
top = gyotoy(gyotoytop)
