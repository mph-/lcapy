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

Note, Lcapy tries to approximate all numbers with integers.  A floating point representation can be found with the `evalf()` method:

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

This is simply the input noise scaled by the amplifier gain :math:`1 + R_2/R_1`.

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
   
