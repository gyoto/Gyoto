/*
    Copyright 2014 Thibaut Paumard

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

// NODISPLAY implies batch mode
if (get_env("GYOTO_CHECK_NODISPLAY")) {
  batch, 1;
  __xytitles=xytitles; __fma=fma; __winkill=winkill; __pli=pli; __plg=plg;
  __pause=pause; __window=window;
  xytitles = fma = winkill = pli = plg = pause = window = noop;
 }

#include "check-helpers.i"
restore, gyoto;

begin_section, "MPI functionalities";

doing, "Checking whether MPI is built-in";
have_mpi=haveMPI();
output, (have_mpi?" yes":" no");

doing, "Calling MPI_Initiliazed";
if (MPI_Initialized()) error, "MPI should not be initialized yet";
done;

doing, "Calling MPI_Init";
if (MPI_Init() && have_mpi) error, "MPI INIT FAILED";
done;

doing, "Calling MPI_Initiliazed again";
inited=MPI_Initialized();
if (have_mpi && !inited) error, "MPI should be initialized by now";
if (!have_mpi && inited) error, "MPI should not be initializable";
done;

doing, "Reading Scenery";
sc = Scenery(GYOTO_EXAMPLES_DIR+"example-complex-astrobj.xml");
done;

doing, "Spawning workers";
sc, mpispawn=4;
done;

doing, "Sending Scenery to the workers";
sc, mpiclone=;
done;

doing, "Ray-tracing with MPI";
data=sc();
done;

doing, "Terminating workers";
sc, mpispawn=0;
done;

doing, "Integrating Scenery without MPI";
sc, nthreads=4;
data2=sc();
done;

doing, "Comparing results";
diff=data-data2;
ind=where(data);
diff(ind)/=data(ind);
mdiff=max(abs(diff));
if (mdiff > 1e-6) error, "Results differ";
output, " OK (max rel. dif.: "+pr1(mdiff)+")";

fov=sc.screen().fov();
npix=sc.screen().resolution();
delta= fov/double(npix);

xx=((indgen(npix)-0.5)*delta-fov/2.)(, -:1:32);
yy=transpose(xx);

sc, mpispawn=8, mpiclone=;

doing, "Integrating whole field, specifying angles, with MPI";
data2 = sc(-xx, yy, );
done;

sc, mpispawn=0;

doing, "Comparing results";
diff=data-data2;
ind=where(data);
diff(ind)/=data(ind);
mdiff=max(abs(diff));
if (mdiff > 1e-6) error, "Results differ";
output, " OK (max rel. dif.: "+pr1(mdiff)+")";

doing, "Deleting Scenery";
sc=[];
done;

doing, "Calling MPI_Finalized";
if (MPI_Finalized()) error, "MPI should not be finalized yet";
done;

doing, "Calling MPI_Finalize";
if (MPI_Finalize() && have_mpi) error, "MPI FINALIZE FAILED";
done;

doing, "Calling MPI_Finalized again";
finited=MPI_Finalized();
if (have_mpi && !finited) error, "MPI should be finalized by now";
if (!have_mpi && finited) error, "MPI should not be finalizable";
done;

doing, "Cleaning";
data=data2=have_mpi=mdiff=diff=inited=finited=[];
done;

end_section, "MPI functionalities";

if (anyof(get_argv() == "check-mpi.i")) quit;
