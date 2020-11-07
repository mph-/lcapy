=========
Tutorials
=========

Opamps
======

An ideal opamp is represented by a voltage controlled voltage source,  The netlist has the form

   >>> E out+ gnd opamp in+ in-  Ac  Ad

Here `Ac` is the open-loop differential gain and `Ad` is the open-loop common-mode gain (zero default).
   

Non-inverting amplifier
-----------------------

   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""
   ... E 1 0 opamp 3 2 Ad; right
   ... W 2 2_1; down
   ... R1 2_1 0_2 R1; down
   ... R2 2_1 1_1 R2; right
   ... W 1 1_1; down
   ... W 3_1 3_2; down
   ... Vs 3_2 0_3 Vs; down
   ... W 0_3 0_1; down
   ... W 3_1 3_3; right
   ... W 3_3 3; right
   ... W 1 1_2; right
   ... P 1_2 0; down
   ... W 0_1 0_2; right
   ... W 0_2 0; right""")
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary
   >>> a.draw()

.. image:: examples/tutorials/opamps/opamp-noninverting-amplifier1.png
   :width: 12cm

The output voltage (at node 1) is found using::           

   >>> Vo = a[1].V(t)
   >>> Vo.pprint()
   Vₛ⋅(A_d⋅R₁ + A_d⋅R₂)
   ────────────────────
      A_d⋅R₁ + R₁ + R₂

When the open-loop differential gain is infinite, the gain just depends on the resistor values::      
           
   >>> Vo.limit('Ad', oo).pprint()
   Vₛ⋅(R₁ + R₂)
   ────────────
        R₁     

Let's now add some common-mode gain to the opamp by overriding the `E` component::

   >>> b = a.copy()
   >>> b.add('E 1 0 opamp 3 2 Ad Ac; right')
           
The output voltage (at node 1) is now::           

   >>> Vo = b[1].V(t)
   >>> Vo.pprint()
   Vₛ⋅(A_c⋅R₁ + A_c⋅R₂ + 2⋅A_d⋅R₁ + 2⋅A_d⋅R₂)
   ──────────────────────────────────────────
        -A_c⋅R₁ + 2⋅A_d⋅R₁ + 2⋅R₁ + 2⋅R₂

Setting the open-loop common-mode gain to zero gives the same result as before::
        
   >>> Vo.limit('Ac', 0).pprint()
   A_d⋅R₁⋅Vₛ + A_d⋅R₂⋅Vₛ
   ─────────────────────
      A_d⋅R₁ + R₁ + R₂

When the open-loop differential gain is infinite, the common-mode gain
of the opamp is insignificant and the gain of the amplifier just
depends on the resistor values::
           
   >>> Vo.limit('Ad', oo).pprint()
   Vₛ⋅(R₁ + R₂)
   ────────────
        R₁     

Let's now consider the input impedance of the amplifier::

   >>> a.impedance(3, 0).pprint()
   0

This is not the desired answer since node 3 is shorted to node 0 by the voltage source `Vs`.  If we try to remove this, we get::

   >>> c = a.copy()
   >>> c.remove('Vs')
   >>> c.impedance(3, 0).pprint()  
   ValueError: The MNA A matrix is not invertible for time analysis because:
   1. there may be capacitors in series;
   2. a voltage source might be short-circuited;
   3. a current source might be open-circuited;
   4. a dc current source is connected to a capacitor (use step current source).
   5. part of the circuit is not referenced to ground

In this case it is reason 3.  This is because Lcapy connects a 1\,A current source across nodes 3 and 0 and tries to measure the voltage to determine the impedance.   However, node 3 is floating since an ideal opamp has infinite input impedance.  To keep Lcapy happy, we can explicitly add a resistor between nodes 3 and 0,

   >>> c.add('Rin 3 0')
   >>> c.impedance(3, 0).pprint()

Now, not surprisingly,

   >>> c.impedance(3, 0).pprint()  
   Rᵢₙ


Inverting amplifier
-------------------

   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""   
   ... E 1 0 opamp 3 2 Ad; right, flipud
   ... R1 4 2; right
   ... R2 2_2 1_1; right
   ... W 2 2_2; up
   ... W 1 1_1; up
   ... W 4 4_2; down=0.5
   ... Vs 4_2 0_3; down
   ... W 0_3 0_1; down=0.5
   ... W 3 0_2; down
   ... W 1 1_2; right
   ... P 1_2 0; down
   ... W 0_1 0_2; right
   ... W 0_2 0; right
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary""")
   >>> a.draw()

.. image:: examples/tutorials/opamps/opamp-inverting-amplifier1.png
   :width: 12cm

The output voltage (at node 1) is found using::

  >>> Vo = a[1].V(t)
  >>> Vo.pprint()
    -A_d⋅R₂⋅vₛ(t)  
   ────────────────
   A_d⋅R₁ + R₁ + R₂

In the limit when the open-loop differential gain is infinite gain of
the amplifier just depends on the resistor values::
           
   >>> Vo.limit('Ad', oo).pprint()
   -R₂⋅vₛ(t) 
   ──────────
       R₁    

Note, the output voltage is inverted compared to the source voltage.

The input impedance can be found by removing `Vs`::

   >>> a.remove('Vs')
   >>> a.impedance(4, 0)
   A_d⋅R₁ + R₁ + R₂
   ────────────────
       A_d + 1     

In the limit with infinite open-loop differential gain::

   >>> a.impedance(4,0).limit('Ad', oo)                                            R₁

   

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

Note, this is the (one-sided) amplitude spectral density with units of volts per root hertz.  Here `T` is the absolute temperature in degrees kelvin, `k` is Boltzmann's constant, and :math:`\omega` is the angular frequency.  The expression can be made a function of linear frequency using:

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


Opamp non-inverting amplifier
-----------------------------

This tutorial looks at the noise from an opamp non-inverting
amplifier.  It uses an ideal opamp with open-loop gain `A` augmented
with a voltage source representing the input-referred opamp voltage
noise, and current sources representing the input-referred opamp
current noise.

   >>> from lcapy import *
   >>> a = Circuit("""
   ... Rs 1 0; down
   ... Vn 1 2 noise; right
   ... W 2 3; right
   ... In1 2 0_2 noise; down, l=I_{n+}
   ... W 0 0_2; right
   ... In2 5 0_5 noise; down, l=I_{n-}
   ... W 5 4; right
   ... W 0_2 0_5; right
   ... W 4 6; down
   ... R1 6 0_6; down
   ... W 0_5 0_6; right
   ... R2 6 7; right
   ... W 8 7; down
   ... E 8 0 opamp 3 4 A; right
   ... W 8 9; right
   ... W 0_6 0_9; right
   ... P 9 0_9; down
   ... ; draw_nodes=connections, label_nodes=none""")
   >>> a.draw()

.. image:: examples/tutorials/opampnoise/opamp-noninverting-amplifier.png
   :width: 10cm

The noise ASD at the input of the opamp is
           
   >>> a[3].V.n
      ____________________________
     ╱    ⎛   2           ⎞     2 
   ╲╱  Rₛ⋅⎝Iₙ₁ ⋅Rₛ + 4⋅T⋅k⎠ + Vₙ  

This is independent of frequency and thus is white.  In practice, the voltage and current noise of an opamp has a 1/f component at low frequencies.

The noise at the output of the amplifier is

   >>> a[8].V.n   
        _____________________________________________________
       ╱    2   2          2      2   2   2     2          2 
   A⋅╲╱  Iₙ₁ ⋅Rₛ ⋅(R₁ + R₂)  + Iₙ₂ ⋅R₁ ⋅R₂  + Vₙ ⋅(R₁ + R₂)  
   ──────────────────────────────────────────────────────────
                         A⋅R₁ + R₁ + R₂                      

Assuming an infinite open-loop gain this simplifies to

   >>> a[8].V.n.limit('A', oo)
      _____________________________________________________
     ╱    2   2          2      2   2   2     2          2 
   ╲╱  Iₙ₁ ⋅Rₛ ⋅(R₁ + R₂)  + Iₙ₂ ⋅R₁ ⋅R₂  + Vₙ ⋅(R₁ + R₂)  
   ────────────────────────────────────────────────────────
                              R₁                           

This is simply the input noise scaled by the amplfier gain :math:`1 + R_2/R_1`.

So far the analysis has ignored the noise due to the feedback resistors.   The noise from these resistors can be modelled with the `noisy()` method of the circuit object.

   >>> b = a.noisy()
   >>> b.draw()

.. image:: examples/tutorials/opampnoise/opamp-noninverting-amplifier-noisy.png
   :width: 10cm


Let's choose :math:`R2 = (G - 1) R_1` where :math:`G` is the closed-loop gain:

   >>> c = b.subs({'R2':'(G - 1) * R1'})
   >>> c[8].V.n.limit('A', oo)

Unfortunately, this becomes unmanageable since SymPy has to assume that :math:`G` may be less than one.   So instead, let's choose :math:`G=10`,

   >>> c = b.subs({'R2':'(10 - 1) * R1'})
   >>> c[8].V.n.limit('A', oo)
      ________________________________________________________________
     ╱        2   2         2   2                                   2 
   ╲╱  100⋅Iₙ₁ ⋅Rₛ  + 81⋅Iₙ₂ ⋅R₁  + 360⋅R₁⋅T⋅k + 400⋅Rₛ⋅T⋅k + 100⋅Vₙ  

In practice, both noise current sources have the same ASD.  Thus

   >>> c = b.subs({'R2':'(10 - 1) * R1', 'In2':'In1'})
   >>> c[8].V.n.limit('A', oo)
      _________________________________________________________________
     ╱        2   2        ⎛     2            ⎞                      2 
   ╲╱  100⋅Iₙ₁ ⋅Rₛ  + 9⋅R₁⋅⎝9⋅Iₙ₁ ⋅R₁ + 40⋅T⋅k⎠ + 400⋅Rₛ⋅T⋅k + 100⋅Vₙ  

The noise is minimised by keeping `R1` as small as possible.  However, for high gains, the noise is dominated by the opamp noise.  Ideally, `Rs` needs to be minimised.  However, if it is large, it is imperative to choose a CMOS opamp with a low noise current.   Unfortunately, these amplifiers have a higher noise voltage than bipolar opamps.
   


Shield guard
============

Electrostatic shields are important to avoid capacitive coupling of intefererence into signals.  However, the capacitance between the signal and cable shields lowers the input impedance of an amplifier.


   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""
   ... Vs 14 12 ac; down
   ... Rs 14 1; right
   ... Cable1; right=4, dashed, kind=coax, l=
   ... W 1 Cable1.in; right=0.5
   ... W Cable1.out 2; right=0.5
   ... W Cable1.ognd 10; down=0.5
   ... Cc Cable1.mid Cable1.b; down=0.2, dashed, scale=0.6
   ... W 2 11; right=0.75
   ... W 7 10; up=0.5
   ... E2 15 0 opamp 11 17 A_1; right, scale=0.5
   ... W 17 18; down
   ... W 12 7; right
   ... W 7 18; right
   ... W 18 0; down=0.2, sground
   ... Rin 11 17; down
   ... ; label_nodes=none, draw_nodes=connections
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary
   >>> a.draw()

.. image:: examples/tutorials/shield-guard/shield-ground.png
   :width: 30cm


To find the input impedance it is first necessary to disconnect the
source, for example,

    >>> a.remove(('Vs', 'Rs'))

The impedance seen across `Rin` can be then found using:

    >>> Z = a.impedance('Rin')           
        1        
    ─────────────────
        ⎛       1   ⎞
    C_c⋅⎜s + ───────⎟
        ⎝    C_c⋅Rᵢₙ⎠

This impedance is the parallel comination of the input resistance Rin and the impedance of the cable capacitance.   Thus at high frequencies the impedance drops.
        

Shield guard circuits are used to mitigate the capacitance between a cable signal and the cable shield.  For example:


   >>> from lcapy import Circuit, t, oo
   >>> b = Circuit("""
   ... Vs 14 12 ac; down
   ... Rs 14 1; right
   ... Cable1; right=4, dashed, kind=coax, l=
   ... W 1 Cable1.in; right=0.5
   ... W Cable1.out 2; right=0.5
   ... W Cable1.ognd 10; down=0.5
   ... Cc Cable1.mid Cable1.b; down=0.2, dashed, scale=0.6
   ... W 2 3; right=1.5
   ... W 3 11; right=0.75
   ... W 3 4; down=0.5
   ... W 5 6; down=0.5
   ... W 6 7; left
   ... W 7 10; up=0.5
   ... R 10 8; right
   ... E1 8 0 opamp 4 5 A_2; left=0.5, mirror, scale=0.5
   ... E2 15 0 opamp 11 17 A_1; right, scale=0.5
   ... W 17 18; down
   ... W 12 18; right
   ... W 18 0; down=0.2, sground
   ... Rin 11 17; down
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary
   >>> a.draw()

.. image:: examples/tutorials/shield-guard/shield-guard.png
   :width: 25cm

Again, to find the input impedance it is first necessary to disconnect the
source, for example,

    >>> b.remove(('Vs', 'Rs'))

The impedance seen across `Rin` can be then found using:

    >>> Z = b.impedance('Rin')
    ⎛Rᵢₙ⋅(A₂ + C_c⋅R⋅s + 1)⎞
    ⎜──────────────────────⎟
    ⎝    C_c⋅(R + Rᵢₙ)     ⎠
    ────────────────────────
               A₂ + 1       
      s + ───────────────   
          C_c⋅R + C_c⋅Rᵢₙ   

However, when the open-loop gain, A2, of the shield-guard amplifier is large then

    >>> Z.limit('A_2', oo)
     Rᵢₙ

Thus the input impedance does not depend on Cc.  In practice, the open-loop gain is not infinite and reduces with frequency and so the guarding does not help at very high frequencies.
