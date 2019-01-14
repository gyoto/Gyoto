import numpy as np
import matplotlib.pyplot as plt
import pyfits as fits
import cv2
import math
import matplotlib.patches as mpatches


image_file="gyoto_image.fits"

fov=250 # field of view in µas, as given in the XML

convertSItoCGS=1e3 # SI to cgs conversion for I_nu (gyoto speaks in SI)

mycmap="YlGnBu" # python cmap choice


###################################################
############ End of input parameters ##############
############  Nothing to edit below  ##############
###################################################

liminf=-fov/2.
limsup=-liminf

image = fits.getdata(image_file)[0]
image*=convertSItoCGS

print("Min,max intensity= ",image.min(),image.max())

plt.figure(0)
plt.clf()
plt.ion()
plt.imshow(image,interpolation='nearest',origin='lower',extent=(-fov/2.,fov/2.,-fov/2.,fov/2.),cmap=mycmap)
plt.xlim(liminf,limsup)
plt.ylim(liminf,limsup)
plt.colorbar()
plt.xlabel("x ($\mu$as)", size=12)
plt.ylabel("y ($\mu$as)", size=12)

NN=image.shape[0]

image_moments=np.zeros((NN, NN, 1), dtype = "float")
image_moments[:,:,0]=image # for some reason, cv2 needs this translation,
# I guess due to some dimension convention
    
allmoms=cv2.moments(image_moments)
mu00=allmoms["m00"]
mu11=allmoms["mu11"]
mu20=allmoms["mu20"]
mu02=allmoms["mu02"]
mup02=mu02/mu00
mup20=mu20/mu00
mup11=mu11/mu00
# ellipse axes in pixels, convert to muas
pxlinmuas=fov/NN
semimajor=math.sqrt(2.*(mup02+mup20
                        +math.sqrt(4.*mup11**2+(mup20-mup02)**2)))
semiminor=math.sqrt(2.*(mup02+mup20
                        -math.sqrt(4.*mup11**2+(mup20-mup02)**2)))
semimajor*=pxlinmuas
semiminor*=pxlinmuas
print("Ellipse typical diameter (muas)= ",2.*0.5*(semimajor+semiminor))
theta=0.5*math.atan2(2.*mup11,(mup20-mup02))*180./math.pi
ellcenxmax=np.where(image==image.max())[1][0]
ellcenymax=np.where(image==image.max())[0][0]
ellcenx=(ellcenxmax-NN/2.)*pxlinmuas
ellceny=(ellcenymax-NN/2.)*pxlinmuas
# the *diameters* of the ellipse should be given to mpatches.Ellipse:
ell = mpatches.Ellipse((ellcenx, ellceny), 2*semimajor, 2*semiminor, theta,
                       edgecolor="black", facecolor='none')
plt.gca().add_patch(ell)

plt.show()
