import numpy as np
import sys

# 1D (just for checking)

def mysimps(yy,xx):
    # this mysimps was compared to scipy.integrate.simps
    # and agrees to perfect accuracy
    nn = len(yy)
    ww = weight1D(nn)
    dx = xx[1]-xx[0]
    mysum=0
    for ii in range(nn):
        mysum+=yy[ii]*ww[ii]
    return dx/3.*mysum

def weight1D(nn):
    # nn i the number of x_i where the integrated function will be called;
    # integ = dx/3*(f_1 + 4*f_2 + 2*f_3 + 4*f_4 + 2*f_5 ... + 4*f_nn-1_ + f_nn)
    # with f_i = f(x_i)
    # This function implements the weights of the f_i above.
    if nn%2!=1:
        sys.exit("nn should be odd")
    if nn<5:
        sys.exit("choose at least 5 evaluation points")
    ww = np.zeros(nn)+1.
    for ii in range(1,nn-1): # [1,2,...nn-2]
        if ii%2==1:
            ww[ii]=4
        else:
            ww[ii]=2
    return ww

# 2D

def mysimps2D(image,dx,dy):
    # 2D sampled Simpson integration, does not exist built-in...
    # dx and dy are the (constant) increment x_i+1_ - x_i
    # and y_i+1_ - y_i
    nn = len(image)
    # Simpson needs odd number of samples, so if nn is even
    # check that last col and last line of image are zero
    # (should typically be in Gyoto...) and remove them
    if nn%2==0:
        if np.count_nonzero(image[:,nn-1])==0 \
           and np.count_nonzero(image[nn-1,0])==0:
            image=image[0:nn-1,0:nn-1]
        else:
            sys.exit("The Gyoto resolution (Npix) should be odd, \
            or the last row and column should be zero")
    nn = len(image)
    ww = weight2D(nn)
    mysum=0
    for ii in range(nn):
        for jj in range(nn):
            mysum+=image[ii,jj]*ww[ii,jj]
    #mysum = sum(map(sum,image*ww))
    return dx*dy/9.*mysum    

def weight2D(nn):
    # Produce the weight matrix:
    # 1   4   2   4  ... 2   4   1
    # 4   16  8   16 ... 8   16  4
    # 2   8   4   8  ... 4   8   2
    # 4   16  8   16 ... 8   16  4
    # .....
    # 2   8   4   8  ... 4   8   2
    # 4   16  8   16 ... 8   16  4
    # 1   4   2   4  ... 2   4   1
    #
    # Remember that M[ii,jj] in python means line ii, col jj
    if nn%2!=1:
        sys.exit("nn should be odd")
    if nn<5:
        sys.exit("choose at least 5 evaluation points")
    ww = np.zeros((nn,nn))
    for ii in range(nn):
        ww[ii,0]=1
        ww[ii,nn-1]=1
        for jj in range(1,nn-1):
            if jj%2==1:
                ww[ii,jj]=4
            else:
                ww[ii,jj]=2
        if ii%2==1:
            ww[ii,:]*=4
        if ii%2==0 and ii>0 and ii<nn-1:
            ww[ii,:]*=2
    return ww
