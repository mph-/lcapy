=========
Tutorials
=========


Noise analysis
==============


RC network
----------

Let's consider the thermal noise produced by a resistor in parallel
with a capacitor.  Now only the resistor produces thermal noise (white
Gaussian noise) but the capacitor and resistor form a filter so the
resultant noise is no longer white.  The interesting thing is that the
resultant noise voltage only depends on the capacitor value.  This can
be demonstrated using Lcapy.   Let's start by creating the circuit:

   >>> from lcapy import *
   >>> a = Circuit("""
   ... R 1 0; down
   ... W 1 2; right
   ... C 2 0_2; down
   ... W 0 0_2; right""")
   >>> a.draw()

.. image:: examples/tutorials/RCnoise/RCparallel1.png
   :width: 4cm

A noisy circuit model can be created with the `noisy()` method of the circuit object.   For every resistor in the circuit, a noisy voltage source is added in series.  For example,

   >>> b = a.noisy()
   >>> b.draw()

.. image:: examples/tutorials/RCnoise/RCparallel1noisy.png
   :width: 4cm
        
The noise voltage across the capacitor can be found using:

   >>> Vn = b.C.V.n
   >>> Vn
       2⋅√R⋅√T⋅√k   
   ─────────────────
      ______________
     ╱  2  2  2     
   ╲╱  C ⋅R ⋅ω  + 1 

Note, this is the amplitude spectral density with units of volts per root hertz.  Here `T` is the absolute temperature in degrees kelvin, `k` is Boltzmann's constant, and :math:`\omega` is the angular frequency.  The expression can be made a function of linear frequency using:

   >>> Vn(f)
         2⋅√R⋅√T⋅√k      
   ──────────────────────
      ___________________
     ╱    2  2  2  2     
   ╲╱  4⋅π ⋅C ⋅R ⋅f  + 1 

This expression can be plotted if we substitute the symbols with numbers.  Let's choose :math:`T = 293` K, :math:`R = 10` kohm, and :math:`C = 100` nF.

   >>> Vns = Vn.subs({'R':10e3, 'C':100e-9, 'T':293, 'k':1.38e-23})
   >>> Vns(f)
              √101085           
   ─────────────────────────────
                    ____________
                   ╱  2  2      
                  ╱  π ⋅f       
   25000000000⋅  ╱   ────── + 1 
               ╲╱    250000     

Note, Lcapy tries to approximate all numbers with integers.  A floating point representation can be found with the `evalf()` method:

   >>> Vns(f).evalf()               
                                                     -0.5
                       ⎛                     2      ⎞    
   1.27175469332729e-8⋅⎝3.94784176043574e-5⋅f  + 1.0⎠    

The amplitude spectral density of the noise can be plotted by definining a vector of frequency samples:

   >>> from numpy import linspace
   >>> vf = linspace(0, 10e3, 200)
   >>> (Vns(f) * 1e9).plot(vf, plot_type='mag', ylabel='ASD (nV/rootHz'))
 

.. image:: examples/tutorials/RCnoise/RCparallel1noiseplot1.png
   :width: 10cm   

Finally, the rms noise voltage can be found using the `rms()` method.  This integrates the square of the ASD (the power spectral density) over all frequencies and takes the square root.  For this example, the rms value does not depend on R.

   >>> Vn.rms()
   √T⋅√k
   ─────
     √C 
