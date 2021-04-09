=========
Tutorials
=========

Basic circuit theory
====================

When learning circuit theory, the key initial concepts are:

1. Ohm's law
2. Kirchhoff's current law (KCL)
3. Kirchhoff's voltage law (KVL)
4. Superposition
5. Norton and Thevenin transformations


DC voltage divider
------------------
   
Consider the DC voltage divider circuit defined by::
  
    >>> from lcapy import Circuit
    >>> a = Circuit("""
    ... V 1 0 6; down=1.5
    ... R1 1 2 2; right=1.5
    ... R2 2 0_2 4; down
    ... W 0 0_2; right""")
    >>> a.draw()
                    
.. image:: examples/tutorials/basic/VRR1.png
   :width: 6cm   

The voltage at node 1 (with respect to the ground node 0) is defined by the voltage source::

   >>> a.V.V
   6

The total resistance is::

   >>> a.R1.R + a.R2.R
   6

and thus using Ohm's law the current through R1 and R2 is::

   >>> I = a.V.V / (a.R1.R + a.R2.R)
   >>> I
   1

Thus again using Ohm's law, the voltage across R2 is::
  
   >>> I * a1.R2.R
   4

Of course, these values can be automatically calculated using Lcapy.  For example,
the voltage at node 2 (with respect to the ground node 0) is::

   >>> a[2].V
   4

This is equivalent to the voltage across R2::
           
   >>> a.R2.V
   4   

The current through R1 is::

   >>> a.R1.I
   1

From Kirchhoff's current law, this is equivalent to the current though R2 and V::

   >>> a.R2.I
   1

   >>> a.V.I
   1     
   
The general result can be obtained by evaluating this circuit symbolically::

    >>> from lcapy import Circuit
    >>> a = Circuit("""
    ... V 1 0 dc
    ... R1 1 2
    ... R2 2 0_2
    ... W 0 0_2""")
    >>> a.R2.V
     R₂⋅V   
    ───────
    R₁ + R₂  

Note, the keyword dc is required here for the voltage source otherwise an arbitrary voltage
source is assumed.
     

AC (phasor) analysis of RC circuit
----------------------------------
   
Consider the circuit defined by::
  
    >>> from lcapy import Circuit
    >>> a = Circuit("""
    ... V 1 0 ac 6; down=1.5
    ... R 1 2 2; right=1.5
    ... C 2 0_2 4; down
    ... W 0 0_2; right""")
    >>> a.draw()
                    
.. image:: examples/tutorials/basic/VRC1.png
   :width: 6cm   

Here the ac keyword specifies that the voltage source is a phasor of angular frequency :math:`\omega_0`.

The voltage across the voltage source is given using::

    >>> a.V.V
    {ω₀: 6}

This indicates a phasor of angular frequency :math:`\omega_0` with an amplitude 6 V.  The time domain representation is::

    >>> a.V.V(t)  
    6⋅cos(ω₀⋅t)

The voltage across the capacitor is::

    >>> a.C.V
   ⎧    ⎛-3⋅ⅉ ⎞⎫
   ⎪    ⎜─────⎟⎪
   ⎪    ⎝  4  ⎠⎪
   ⎨ω₀: ───────⎬
   ⎪          ⅉ⎪
   ⎪     ω₀ - ─⎪
   ⎩          8⎭
 

This can be simplified::

    >>> a.C.V.simplify()
   ⎧     -6⋅ⅉ   ⎫
   ⎨ω₀: ────────⎬
   ⎩    8⋅ω₀ - ⅉ⎭


The magnitude of the phasor voltage is::

    >>> a.C.V.magnitude
       ___________________
      ╱          2        
    ╲╱  589824⋅ω₀  + 9216 
    ──────────────────────
               2          
        1024⋅ω₀  + 16     

and the phase is::

    >>> a.C.V.phase
    -atan(8⋅ω₀)
  
Finally, the time-domain voltage across the capacitor is::

    >>> a.C.V(t)
    48⋅ω₀⋅sin(ω₀⋅t)   6⋅cos(ω₀⋅t)
    ─────────────── + ───────────
           2               2    
       64⋅ω₀  + 1      64⋅ω₀  + 1


Laplace analysis of RC low-pass filter
--------------------------------------

The following netlist describes a first-order RC low-pass filter (the
`P` components define the input and output ports)::

    >>> from lcapy import Circuit
    >>> a = Circuit("""
    ... P1 1 0; down=1.5, v=v_i(t)
    ... R 1 2 2; right=1.5
    ... C 2 0_2 {1/4}; down
    ... W 0 0_2; right
    ... W 2 3; right
    ... W 0_2 0_3; right
    ... P2 3 0_3; down, v^=v_o(t)"""
    >>> a.draw()
                    
.. image:: examples/tutorials/basic/VRC2.png
   :width: 6cm

Here :math:`v_i(t)` is the input voltage and :math:`v_o(t)` is the output voltage.  The transfer function of the filter can be found by specifying nodes:

   >>> H = a.transfer(1, 0, 3, 0)

or components:

   >>> H = a.P1.transfer('P2')

In both cases, the transfer function is::

   >>> H
      s     
   ───────
   (s + 2)

For the input signal, let's consider a sinewave of angular frequency 3 rad/s that switches 'on' at :math:`t=0`::

   >>> v_i = voltage(sin(3 * t) * u(t))

The output voltage can be found by connecting a voltage source with
this signal to the circuit and using Lcapy to find the result.
However, let's use Laplace transforms to find the result.  For this signal, its Laplace transform is::

   >>> V_i = v_i(s)
   >>> V_i
      3   
    ──────
     2    
    s  + 9

The Laplace transform of the output voltage is found by multiplying this with the transfer function::

   >>> V_o = V_i * H
   >>> V_o
             6          
   ────────────────────
    3      2           
   s  + 2⋅s  + 9⋅s + 18

This has three poles: two from the input signal and one from the transfer function of the filter.  This can be seen from the zero-pole-gain form of the response:

   >>> V_o.ZPK()
                   6             
   ───────────────────────────
   (s + 2)⋅(s - 3⋅ⅉ)⋅(s + 3⋅ⅉ)


Using an inverse Laplace transform, the output voltage signal in the time-domain is::

   >>> v_o = V_o(t)
   >>> v_o
     ⎛                                               -2⋅t⎞     
     ⎜2⋅sin(3⋅t)   cos(3⋅t)   (-2 - 3⋅ⅉ)⋅(-2 + 3⋅ⅉ)⋅ℯ    ⎟     
   6⋅⎜────────── - ──────── + ─────────────────────────── ⎟⋅u(t)
     ⎝    39          13                  169            ⎠     

This can be simplified, however, SymPy has trouble with this as a whole.  Instead it is better
to simplify the expression term by term::

  >>> v_o.simplify_terms()
                                          -2⋅t     
   4⋅sin(3⋅t)⋅u(t)   6⋅cos(3⋅t)⋅u(t)   6⋅ℯ    ⋅u(t)
   ─────────────── - ─────────────── + ────────────
          13                13              13     

The first two terms represent the steady-state reponse and the third
term represents the transient response due to the sinewave switching
'on' at :math:`t=0`.  The steady-state response is the sum of a
sinewave and cosinewave of the same frequency; this is equivalent to a
phase-shifted sinewave.  This can be seen using the `simplify_sin_cos`
method::

   >>> v_o.simplify_sin_cos(as_sin=True)

            ⎛      π            ⎞                    
   2⋅√13⋅sin⎜3⋅t - ─ + atan(2/3)⎟⋅u(t)      -2⋅t     
            ⎝      2            ⎠        6⋅ℯ    ⋅u(t)
   ─────────────────────────────────── + ────────────
                    13                        13     

Here the phase delay is `-pi/2 + atan(2/3)` or about -56 degrees::

  >>> ((-pi/2 + atan(2/3)) / pi * 180).fval
  -56.3

The input and output signals can be plotted using::

   >>> ax = v_i.plot((-1, 10), label='input')
   >>> ax = v_o.plot((-1, 10), axes=as, label='output')
   >>> ax.legend()

.. image:: examples/tutorials/basic/VRC2plot.png
   :width: 12cm

Notice the effect of the transient at the start before the response tends to the steady state response.

The phase response of the filter can be plotted as follows::

  >>> H(jw).phase_degrees.plot((0, 10))

.. image:: examples/tutorials/basic/VRC2Hphaseplot.png
   :width: 12cm

Notice that the phase shift is -56 degrees at an angular frequency of 3 rad/s.   

The amplitude response of the filter can be plotted as follows::

  >>> H(jw).magnitude.plot((0, 10))

.. image:: examples/tutorials/basic/VRC2Hmagplot.png
   :width: 12cm

For a Bode plot, the angular frequency is plotted on a logarithmic
scale and the amplitude response is plotted in dB::

  >>> H(jw).dB.plot((0, 10), log_frequency=True)

.. image:: examples/tutorials/basic/VRC2HdBplot.png
   :width: 12cm                        


           

Superposition of AC and DC
--------------------------

Here's an example circuit comprised of two AC sources and a DC source:
    >>> from lcapy import Circuit
    >>> a = Circuit("""
    ... V1 1 0 {2 * sin(3*t)}; down=1.5
    ... V2 1 2 {3 * cos(4*t)}; right=1.5
    ... V3 3 2 4; left=1.5
    ... R 3 0_3; down
    ... W 0 0_3; right""")
    >>> a.draw()
                    
.. image:: examples/tutorials/basic/Vsup3.png
   :width: 9cm   

The voltage across the resistor is the sum of the three voltage
sources.  This is shown as a superposition::

    >>> a.R.V
    {dc: 4, 3: -2⋅ⅉ, 4: -3}

This shows that there is a DC component of 4 V added to two phasors;
one with an angular frequency of 3 rad/s and the other with angular
frequency of 4 rad/s.

There are a number of ways in which the signal components can be extracted.
For example, the phase of the 3 rad/s phasor can be found using::

    >>> a.R.V[3].phase
    -π 
    ───
     2

Similarly, the magnitude of of the 4 rad/s phasor can be found using::

    >>> a.R.V[4].magnitude
    3

The DC component can be extracted using::

    >>> a.R.V.dc
    4    
  
Alternatively, since DC is a phasor of angular frequency 0 rad/s::

   >>> a.R.V[0]
   4  
  
The overall time varying voltage can be found using::

   >>> a.R.V(t)
   2⋅sin(3⋅t) - 3⋅cos(4⋅t) + 4

    
       
Initial value problem
=====================

Consider the series R-L-C circuit described by the netlist:

.. literalinclude:: examples/tutorials/ivp/circuit-RLC-ivp1.sch

Note, to specify the initial conditions, the capacitance and
inductance values must be explicitly defined.
                    
This can be loaded by Lcapy and drawn using:

    >>> from lcapy import Circuit, s, t
    >>> a = Circuit("circuit-RLC-ivp1.sch")
    >>> a.draw()
                   
 

This circuit has a specified initial voltage for the capacitor and a
specified initial current for the inductor.  Thus, it is solved as an
initial value problem.  This will give the transient response for
:math:`t \ge 0`.  Note, the initial values usually arise for switching
problems where the circuit topology changes. 

The s-domain voltage across the resistor can be found using::

   >>> c.R.V(s)
   ⎛R⋅(L⋅i₀⋅s - v₀)⎞
   ⎜───────────────⎟
   ⎝       L       ⎠
   ─────────────────
      2   R⋅s    1  
     s  + ─── + ─── 
           L    C⋅L 

This can be split into terms, one for each initial value, using::

   >>> a.R.V(s).expandcanonical() 
        R⋅i₀⋅s              R⋅v₀       
   ────────────── - ──────────────────
    2   R⋅s    1      ⎛ 2   R⋅s    1 ⎞
   s  + ─── + ───   L⋅⎜s  + ─── + ───⎟
         L    C⋅L     ⎝      L    C⋅L⎠

Lcapy can convert this expression into the time-domain but the result
is complicated.  This is because SymPy does not know how to simplify
the expression since it cannot tell if the poles are complex
conjugates, distinct real, or repeated real.  Let's have a look at the
poles::

  >>> a.R.V(s).poles()
  ⎧           ____________                ____________   ⎫
  ⎪          ╱    2                      ╱    2          ⎪
  ⎨   R    ╲╱  C⋅R  - 4⋅L         R    ╲╱  C⋅R  - 4⋅L    ⎬
  ⎪- ─── + ───────────────: 1, - ─── - ───────────────: 1⎪
  ⎩  2⋅L        2⋅√C⋅L           2⋅L        2⋅√C⋅L       ⎭

Thus it can be seen that if :math:`C R^2 \ge 4 L` then the poles are
real otherwise they are complex.

To get a simpler result that does not depend on the unknown component
values, let's parameterize the expression for the voltage across the
resistor::

   >>> VR, defs = a.R.V(s).parameterize()
   >>> VR
         K⋅(L⋅i₀⋅s - v₀)      
   ──────────────────────────
        ⎛  2               2⎞
   L⋅i₀⋅⎝ω₀  + 2⋅ω₀⋅s⋅ζ + s ⎠
   >>> defs
   ⎧                    1          √C⋅R⎫
   ⎨K: R⋅i₀, omega_0: ─────, zeta: ────⎬
   ⎩                  √C⋅√L        2⋅√L⎭

Unfortunately, converting VR into the time-domain also results in a
complicated expression that SymPy cannot simplify.  Instead, it is
better to use an alternative parameterization::

   >>> VR, defs = a.R.V(s).parameterize(zeta=False)
   >>> VR
          K⋅(L⋅i₀⋅s - v₀)        
   ──────────────────────────────
        ⎛  2    2              2⎞
   L⋅i₀⋅⎝ω₁  + s  + 2⋅s⋅σ₁ + σ₁ ⎠

   >>> defs
   ⎧                     ______________              ⎫
   ⎪                    ╱      2                     ⎪
   ⎨                  ╲╱  - C⋅R  + 4⋅L             R ⎬
   ⎪K: R⋅i₀, omega_1: ─────────────────, sigma_1: ───⎪
   ⎩                        2⋅√C⋅L                2⋅L⎭


The time-domain response can now be found::

   >>> VR(t)
     ⎛     ⎛                       -σ₁⋅t          ⎞       -σ₁⋅t          ⎞           
     ⎜     ⎜ -σ₁⋅t             σ₁⋅ℯ     ⋅sin(ω₁⋅t)⎟   v₀⋅ℯ     ⋅sin(ω₁⋅t)⎟           
   K⋅⎜L⋅i₀⋅⎜ℯ     ⋅cos(ω₁⋅t) - ───────────────────⎟ - ───────────────────⎟           
     ⎝     ⎝                            ω₁        ⎠            ω₁        ⎠           
   ───────────────────────────────────────────────────────────────────────  for t ≥ 0
                                     L⋅i₀                                            
                                                                                     

Finally, the result in terms of R, L, and C can be found by
substituting the parameter definitions::

   >>> VR(t).subs(defs)
   
However, the result is too long to show.   
   

Opamps
======

An ideal opamp is represented by a voltage controlled voltage source,  The netlist has the form

   >>> E out+ gnd opamp in+ in-  Ad  Ac

Here `Ad` is the open-loop differential gain and `Ac` is the open-loop common-mode gain (zero default).
   

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

In this case it is reason 3.  This is because Lcapy connects a 1 A current source across nodes 3 and 0 and tries to measure the voltage to determine the impedance.   However, node 3 is floating since an ideal opamp has infinite input impedance.  To keep Lcapy happy, we can explicitly add a resistor between nodes 3 and 0,

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

In the limit when the open-loop differential gain is infinite, the gain of
the amplifier just depends on the resistor values::
           
   >>> Vo.limit('Ad', oo).pprint()
   -R₂⋅vₛ(t) 
   ──────────
       R₁    

Note, the output voltage is inverted compared to the source voltage.

The input impedance can be found by removing `Vs` (since it has zero resistance)::

   >>> a.remove('Vs')
   >>> a.impedance(4, 0)
   A_d⋅R₁ + R₁ + R₂
   ────────────────
       A_d + 1     

In the limit with infinite open-loop differential gain::

   >>> a.impedance(4,0).limit('Ad', oo)                                            R₁


However, in practice, the open-loop gain decreases with frequency and so at high frequencies,

   >>> a.impedance(4,0).limit('Ad', 0)
   R₁ + R₂
   

Transimpedance amplifier
------------------------

   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""   
   ...E 1 0 opamp 3 2 Ad; right, flipud
   ...W 4 2; right
   ...R 2_2 1_1; right
   ...W 2 2_2; up
   ...W 1 1_1; up
   ...W 4 4_2; down=0.5
   ...Is 4_2 0_3; down
   ...W 0_3 0_1; down=0.5
   ...W 3 0_2; down
   ...W 1 1_2; right
   ...P 1_2 0; down
   ...W 0_1 0_2; right
   ...W 0_2 0; right
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary""")
   >>> a.draw()

.. image:: examples/tutorials/opamps/opamp-transimpedance-amplifier1.png
   :width: 12cm

The output voltage (at node 1) is found using::

  >>> Vo = a[1].V(t)
  >>> Vo.pprint()
  -A_d⋅R⋅iₛ(t) 
  ─────────────
     A_d + 1   

In the limit when the open-loop differential gain is infinite gain, the gain of
the amplifier just depends on the resistor value::
           
   >>> Vo.limit('Ad', oo).pprint()
   -R⋅iₛ(t) 

Note, the output voltage is inverted compared to the source current.

The input impedance can be found using (the ideal current source has infinite resistance and does not need removing)::

   >>> a.impedance(4, 0)
      R   
   ───────
   A_d + 1

In the limit with infinite open-loop differential gain::

   >>> a.impedance(4,0).limit('Ad', oo)                                            0

However, in practice, the open-loop gain decreases with frequency and so at high frequencies,

   >>> a.impedance(4,0).limit('Ad', 0)
   R

Piezo transducer amplifier
--------------------------

Consider the opamp noise model for a piezo transducer amplifier (note,
the thermal noise of the resistors is ignored):

   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""   
   ... Cs 1 0; down=4
   ... W 1 1_1; right
   ... Rs 1_1 0_1; down=4
   ... W 0 0_1; right
   ... W 1_1 1_2; right=2
   ... Vn 1_2 2 noise; right
   ... E 5_1 0 opamp 2 3_2 A; right
   ... W 3 3_1; right
   ... W 3_1 3_2; right
   ... W 3 3_3; down
   ... R1 3_3 4; down
   ... C 4 0_2; down
   ... W 0_1 0_2; right
   ... W 3_1 3_4; down=1.5
   ... Inn 3_4 0_3 noise; down
   ... W 0_2 0_3; right
   ... W 2 2_1; down=2
   ... Inp 2_1 0_4 noise; down
   ... W 0_3 0_4; right
   ... W 3_3 3_5; right=3
   ... R2 3_5 5_2; right
   ... W 5_1 5_2; down=2
   ... W 5_1 5; right
   ... Po 5 0_5; down, v=v_o
   ... W 0_4 0_5; right
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary""")
   >>> a.draw()

.. image:: examples/tutorials/opamps/opamp-piezo-amplifier1.png
   :width: 12cm

This circuit is a non-inverting amplifier where Cs represents the piezo transducer.  The AC voltage gain is set by 1 + R2/R1, for frequencies above the pole created by R1 and C::

   >>> H = a.transfer('Cs', 'Po').limit('A', oo)
   >>> H
   C⋅s⋅(R₁ + R₂) + 1
   ─────────────────
       C⋅R₁⋅s + 1   
   >>> H.ZPK()
             ⎛         1     ⎞
   (R₁ + R₂)⋅⎜s + ───────────⎟
             ⎝    C⋅(R₁ + R₂)⎠
   ───────────────────────────
             ⎛     1  ⎞       
          R₁⋅⎜s + ────⎟       
             ⎝    C⋅R₁⎠       
           
The output noise voltage amplitude spectral density (ASD) can be found using::

  >>> Von = a.Po.V.n(f)

This can be simplfied by assuming an opamp with infinite open-loop gain::

  >>> Von = Von.limit('A', oo)

This expression is still too complicated to print.  This amplifier circuit has zero gain at DC where the noise ASD is::

   >>> Von.limit(f, 0)
      ___________________________
     ╱    2   2      2   2     2 
   ╲╱  Iₙₙ ⋅R₂  + Iₙₚ ⋅Rₛ  + Vₙ

In comparison, the noise ASD at high frequencies is::

   >>> Von.limit(f, oo)
      ________________________________________________
     ╱    2   2   2     2   2             2     2   2 
   ╲╱  Iₙₙ ⋅R₁ ⋅R₂  + R₁ ⋅Vₙ  + 2⋅R₁⋅R₂⋅Vₙ  + R₂ ⋅Vₙ  
   ───────────────────────────────────────────────────
                            R₁                        

Note that this expression does not depend on Rs due to the capacitance, Cs, of the piezo transducer.

Let's choose some values, R1=100 ohm, R2=900 ohm, Cs=1 nF, Rs=100 Mohm, C=100 nF, Vn=2 nV :math:`/\sqrt{\mathrm{Hz}}`, In=5 fA :math:`/\sqrt{\mathrm{Hz}}` (in practice, the current and voltage noise ASD will vary with frequency due to 1/f noise)::

  b = a.subs({'R1':100, 'R2':900, 'Cs':1e-9, 'Rs':100e6, 'C':100e-9, 'Vn':2e-9, 'Inp':5e-15, 'Inn':0.5e-15})

The high frequency noise ASD is::

  >>> b.Po.V.n(f).limit('A', oo).limit(f, oo).evalf(3)
  2.00e-8

This is the opamp noise voltage scaled by the gain of the amplifier (10).  The noise at DC is larger::

  >>> b.Po.V.n(f).limit('A', oo).limit(f, 0).evalf(3)
  5.00e-7  
  
This is due to the current noise of the opamp flowing through Rs.   The noise ASD would be larger if not for C limiting the low-frequency gain of the amplifier.  For example::

   >>> c = b.replace('C', 'W 4 0_2')
   >>> c.Po.V.n(f).limit('A', oo).limit(f, 0).evalf(3)
   5.00e-6

In this case, the DC gain is 10, whereas with C it is 1.

The amplifier output voltage noise ASD and gain can be plotted using::

  >>> from numpy import logspace
  >>> vf = logspace(0, 5, 201)
  >>> Vo = b.Po.V.n(f).limit('A', oo)
  >>> Vo.plot(vf, log_frequency=True, yscale=1e9, ylabel='Voltage noise ASD (nV/rtHz)')  
  >>> H = b.transfer('Cs', 'Po')(f).limit('A', oo)
  >>> H.plot(vf, log_frequency=True, ylabel='Gain (dB)')
  >>> from matplotlib.pyplot import show
  >>> show()

.. image:: examples/tutorials/opamps/opamp-piezo-amplifier2-asd2.png
   :width: 12cm


.. image:: examples/tutorials/opamps/opamp-piezo-amplifier2-frequency-response1.png
   :width: 12cm

So far the thermal noise generated by the resistors has been ignored.
This can be remedied using the `noise_model()` method that adds a noise voltage source in
series with each resistor.  For example,

   >>> r = b.noise_model().subs({'T':293,'k_B':1.38e-23})

Here, a temperature of 293 K (20 C) is assumed and Boltzmann's constant has been substituted.  The high-frequency noise
is increased slightly by the resistors::

   >>> r.Po.V.n(f).limit('A', oo).limit(f, oo).evalf(3)
   2.34e-8


Multiple feedback low-pass filter
---------------------------------

   >>> from lcapy import Circuit, t, oo
   >>> a = Circuit("""
   ... R2 1 2 R2; right
   ... R3 2 3 R3; right
   ... C2 2 0_3 C2; down=1.15
   ... E1 6 0_3 opamp 0 3 A; mirror, scale=0.5, size=0.5, l=A
   ... W 0 0_2; down
   ... W 3 3_1; up=0.75
   ... C3 3_1 6_1 C3; right
   ... W 6_1 6; down
   ... W 2 2_1; up=1
   ... R4 2_1 6_2 R4; right
   ... W 6_2 6_1; down
   ... W 6 7; right
   ... W 0_3 0_2; right
   ... W 0_2 0_7; right
   ... W 0_3 0_1; left=1
   ... P1 1 0_1; down
   ... P2 7 0_7; down
   ... ; draw_nodes=connections, label_ids=none, label_nodes=primary""")
   >>> a.draw()   

.. image:: examples/tutorials/opamps/multiple-feedback-lpf.png
   :width: 12cm

The transfer function of the filter (assuming an ideal opamp with infinite differential open loop gain) can be found using::

  >>> H = a.transfer(1, 0, 7, 0).limit('A', oo)
                                   -R₄                             
   ─────────────────────────────────────────────────────────────
                2                                            
   C₂⋅C₃⋅R₂⋅R₃⋅R₄⋅s  + C₃⋅R₂⋅R₃⋅s + C₃⋅R₂⋅R₄⋅s + C₃⋅R₃⋅R₄⋅s + R₂

The angular frequency response of the filter is::

  >>> Hf = H(jw)
                                    -R₄                                 
   ─────────────────────────────────────────────────────────────────────
                     2                                                  
   - C₂⋅C₃⋅R₂⋅R₃⋅R₄⋅ω  + ⅉ⋅C₃⋅R₂⋅R₃⋅ω + ⅉ⋅C₃⋅R₂⋅R₄⋅ω + ⅉ⋅C₃⋅R₃⋅R₄⋅ω + R₂

   
The DC gain of the filter is::

  >>> Hf(0)
    -R₄ 
   ────
    R₂ 

The second-order response of the transfer function of the filter can
be parameterized::

  >>> Hp, defs = H.parameterize()
  >>> Hp
              K         
   ───────────────────
     2               2
   ω₀  + 2⋅ω₀⋅s⋅ζ + s

where::

   >>> defs['K']
       -1     
   ───────────
   C₂⋅C₃⋅R₂⋅R₃
   >>> defs['omega_0']
                1             
   ───────────────────────────
     ____   ____   ____   ____
   ╲╱ C₂ ⋅╲╱ C₃ ⋅╲╱ R₃ ⋅╲╱ R₄ 

   >>> defs['zeta']
     ____   ____   ____   ____ ⎛  1       1       1  ⎞
   ╲╱ C₂ ⋅╲╱ C₃ ⋅╲╱ R₃ ⋅╲╱ R₄ ⋅⎜───── + ───── + ─────⎟
                               ⎝C₂⋅R₄   C₂⋅R₃   C₂⋅R₂⎠
   ───────────────────────────────────────────────────
                            2                         

The Q of the filter is found using::

  >>> Q = 1 / (defs['zeta'] * 2)
  >>> Q.simplify()
        ____      ____   ____    
      ╲╱ C₂ ⋅R₂⋅╲╱ R₃ ⋅╲╱ R₄     
   ──────────────────────────────
     ____                        
   ╲╱ C₃ ⋅(R₂⋅R₃ + R₂⋅R₄ + R₃⋅R₄)

   
Shield guard
------------

Electrostatic shields are important to avoid capacitive coupling of interference into signals.  However, the capacitance between the signal and cable shields lowers the input impedance of an amplifier.  Consider the circuit described by the following netlist::

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
   :width: 22cm


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

This impedance is the parallel combination of the input resistance Rin and the impedance of the cable capacitance.   Thus at high frequencies the impedance drops.
        

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
    2⋅√R⋅√T⋅╲╱ k_B  
   ─────────────────
      ______________
     ╱  2  2  2     
   ╲╱  C ⋅R ⋅ω  + 1 

Note, this is the (one-sided) amplitude spectral density with units of volts per root hertz.  Here `T` is the absolute temperature in degrees kelvin, `k_B` is Boltzmann's constant, and :math:`\omega` is the angular frequency.  The expression can be made a function of linear frequency using:

   >>> Vn(f)
   2⋅√R⋅√T⋅╲╱ k_B     
   ──────────────────────
      ___________________
     ╱    2  2  2  2     
   ╲╱  4⋅π ⋅C ⋅R ⋅f  + 1 

This expression can be plotted if we substitute the symbols with numbers.  Let's choose :math:`T = 293` K, :math:`R = 10` kohm, and :math:`C = 100` nF.

   >>> Vns = Vn.subs({'R':10e3, 'C':100e-9, 'T':293, 'k_B':1.38e-23})
   >>> Vns(f)
              √101085           
   ─────────────────────────────
                    ____________
                   ╱  2  2      
                  ╱  π ⋅f       
   25000000000⋅  ╱   ────── + 1 
               ╲╱    250000     

Note, Lcapy tries to approximate real numbers with rationals.  A floating point representation can be found with the `evalf()` method:

   >>> Vns(f).evalf()               
                                                     -0.5
                       ⎛                     2      ⎞    
   1.27175469332729e-8⋅⎝3.94784176043574e-5⋅f  + 1.0⎠    

The amplitude spectral density of the noise can be plotted by defining a vector of frequency samples:

   >>> from numpy import linspace
   >>> vf = linspace(0, 10e3, 200)
   >>> (Vns(f) * 1e9).plot(vf, plot_type='mag', ylabel='ASD (nV/rootHz'))
 

.. image:: examples/tutorials/RCnoise/RCparallel1noiseplot1.png
   :width: 10cm   

Finally, the rms noise voltage can be found using the `rms()` method.  This integrates the square of the ASD (the power spectral density) over all frequencies and takes the square root.  For this example, the rms value does not depend on R.

   >>> Vn.rms()
   √T⋅√k_B
   ───────
     √C 

Opamp noise
-----------

The noise introduced by an opamp is characterised by its noise voltage amplitude spectral density (ASD) and it noise current ASD.  At low frequencies the opamp noise is dominated by 1/f noise, also known as flicker noise.   A model for the noise voltage (ASD) is:

.. math::
   \mathcal{V}(f) =  \mathcal{V}_w \sqrt{\frac{f_c}{f} + 1},

where :math:`\mathcal{V}_w` is the noise ASD in the `white` region and
:math:`f_c` is the corner frequency.  For example, a low noise opamp
with a white noise ASD of 1 nV :math:`/\sqrt{\mathrm{Hz}}` and a corner
frequency of 3.5 Hz can be modelled with::

   >>> from lcapy import f, sqrt
   >>> Vn = 1e-9 * sqrt(3.5 / f + 1)
   >>> ax = (Vn * 1e9).plot((-1, 4), loglog=True)
   >>> ax.grid(True, 'both')
   >>> ax.set_ylabel('Noise voltage density nV$/sqrt{\mathrm{Hz}}$')
   >>> ax.set_ylim(0.1, 10)

.. image:: examples/tutorials/opamps/vnoise1.png
   :width: 10cm

The opamp noise current ASD also has flicker noise.  For example::

   >>> from lcapy import f, sqrt
   >>> In = 1e-12 * sqrt(250 / f + 1)
   >>> ax = (In * 1e12).plot((-1, 4), loglog=True)
   >>> ax.grid(True, 'both')
   >>> ax.set_ylabel('Noise voltage density pA$/sqrt{\mathrm{Hz}}$')
   >>> ax.set_ylim(0.1, 10)

.. image:: examples/tutorials/opamps/inoise1.png
   :width: 10cm  


     
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
     ╱    ⎛   2             ⎞    2 
   ╲╱  Rₛ⋅⎝Iₙ₁ ⋅Rₛ + 4⋅T⋅k_B⎠ + Vₙ  

This is independent of frequency and thus is white.  In practice, the voltage and current noise of an opamp has a 1/f component that dominates at low frequencies.

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

This is simply the input noise scaled by the amplifier gain :math:`1 + R_2/R_1`.

So far the analysis has ignored the noise due to the feedback resistors.   The noise from these resistors can be modelled with the `noisy()` method of the circuit object.

   >>> b = a.noisy()
   >>> b.draw()

.. image:: examples/tutorials/opampnoise/opamp-noninverting-amplifier-noisy.png
   :width: 10cm


Let's choose :math:`R2 = (G - 1) R_1` where :math:`G` is the closed-loop gain of the amplifier:

   >>> c = b.subs({'R2':'(G - 1) * R1'})
   >>> c[8].V.n.limit('A', oo)

Unfortunately, this becomes unmanageable since SymPy has to assume that :math:`G` may be less than one.   So instead, let's choose :math:`G=10`,

   >>> c = b.subs({'R2':'(10 - 1) * R1'})
   >>> c[8].V.n.limit('A', oo)
      __________________________________________________________________
     ╱        2   2         2   2                                   2 
   ╲╱  100⋅Iₙ₁ ⋅Rₛ  + 81⋅Iₙ₂ ⋅R₁  + 360⋅R₁⋅T⋅k_B + 400⋅Rₛ⋅T⋅k_B + 100⋅Vₙ  

In practice, both noise current sources have the same ASD.  Thus

   >>> c = b.subs({'R2':'(10 - 1) * R1', 'In2':'In1'})
   >>> c[8].V.n.limit('A', oo)
      _________________________________________________________________
     ╱        2   2        ⎛     2              ⎞                     2 
   ╲╱  100⋅Iₙ₁ ⋅Rₛ  + 9⋅R₁⋅⎝9⋅Iₙ₁ ⋅R₁ + 40⋅T⋅k_B⎠ + 400⋅Rₛ⋅T⋅k + 100⋅Vₙ  

The noise is minimised by keeping `R1` as small as possible.  However, for high gains, the noise is dominated by the opamp noise.  Ideally, `Rs` needs to be minimised.  However, if it is large, it is imperative to choose a CMOS opamp with a low noise current.   Unfortunately, these amplifiers have a higher noise voltage than bipolar opamps.
   
Here's an script that shows the noise contributions due to the opamp
voltage and current noise as well as the thermal noise from the source
and amplifier resistances at 20 degrees C.

    >>> from lcapy import Circuit, sqrt, f, oo
    >>> Rs = 30
    >>> G = 1000
    >>> R1 = 100
    >>> R2 = (G - 1) * R1
    >>> Vn = 1e-9 * sqrt(3.5 / f + 1)
    >>> In = 1e-12 * sqrt(250 / f + 1)
    >>> T = 273 + 20
    >>> k_B = 1.38e-23
    >>> a = Circuit('opamp-noninverting-amplifier.sch')
    >>> an = a.noisy()
    >>> Vno = an[8].V.n(f)
    >>> Vno = Vno.limit('A', oo)
    >>> Vnov = Vno.subs({'R1':R1, 'R2':R2, 'In1':0, 'In2':0, 'Vn':Vn, 'Rs':Rs,
    ... 'k_B':k_B, 'T':0})
    >>> Vnoi = Vno.subs({'R1':R1, 'R2':R2, 'In1':In, 'In2':In, 'Vn':0, 'Rs':Rs,
    ... 'k_B':k_B, 'T':0})
    >>> Vnor = Vno.subs({'R1':R1, 'R2':R2, 'In1':0, 'In2':0, 'Vn':0, 'Rs':Rs,
    ... 'k_B':k_B, 'T':T})
    >>> Vnot = Vno.subs({'R1':R1, 'R2':R2, 'In1':In, 'In2':In, 'Vn':Vn, 'Rs':Rs,
    ... 'k_B':k_B, 'T':T})
    >>> flim = (-1, 4)
    >>> ax = (Vnot * 1e9).plot(flim, loglog=True, label='total')
    >>> ax = (Vnov * 1e9).plot(flim, loglog=True, label='Vn', axes=ax)
    >>> ax = (Vnoi * 1e9).plot(flim, loglog=True, label='In', axes=ax)
    >>> ax = (Vnor * 1e9).plot(flim, loglog=True, label='R', axes=ax)
    >>> ax.set_ylabel('Noise voltage density nV/$\sqrt{\mathrm{Hz}}$')
    >>> ax.grid(True, 'both')
    >>> ax.legend()

.. image:: examples/tutorials/opampnoise/opamp-noninverting-amplifier-noise1.png
   :width: 10cm

In this plot, the blue line denotes the total noise voltage ASD at the
output of the amplifier, orange shows the noise voltage ASD due to the
opamp voltage noise, green shows the noise voltage ASD due to the
opamp noise current flow the source and amplifier resistances, and the
red curve shows the noise voltage ASD due to thermal noise in the
source and amplifier resistances.  In this example, the total noise is
dominated by the thermal noise and thus lower value amplifier
resistances should be selected.
