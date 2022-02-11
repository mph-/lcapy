=========================
Advanced circuit analysis
=========================


Transmission lines
==================

A lossless transmission line with characteristic impedance
:math:`Z_0`, speed of propagation :math:`c`, and length :math:`l`
driven by a step voltage :math:`V_s` through a resistance :math:`R_s`,
and loaded with a resistance :math:`R_l` can be defined as:


.. literalinclude:: examples/tutorials/txline/txline1.sch

Lcapy can draw this using:

   >>> from lcapy import Circuit, s, t
   >>> a = Circuit('txline1.sch')
   >>> a.draw()


.. image:: examples/tutorials/txline/txline1.png
   :width: 8cm

The driving-point impedance looking into the start of the transmission line is:

   >>> a.P1.dpZ
                      ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞⎞
                Rₛ⋅Z₀⋅⎜Rₗ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟⎟
                      ⎝       ⎝ c ⎠          ⎝ c ⎠⎠
   ───────────────────────────────────────────────────────────────────
             ⎛l⋅s⎞             ⎛l⋅s⎞             ⎛l⋅s⎞     2     ⎛l⋅s⎞
   Rₗ⋅Rₛ⋅sinh⎜───⎟ + Rₗ⋅Z₀⋅cosh⎜───⎟ + Rₛ⋅Z₀⋅cosh⎜───⎟ + Z₀ ⋅sinh⎜───⎟
             ⎝ c ⎠             ⎝ c ⎠             ⎝ c ⎠           ⎝ c ⎠

Similarly, the driving-point impedance looking into the end of the transmission line is:

   >>> a.P2.dpZ
                         ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞⎞
                   Rₗ⋅Z₀⋅⎜Rₛ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟⎟
                         ⎝       ⎝ c ⎠          ⎝ c ⎠⎠
   ───────────────────────────────────────────────────────────────────
             ⎛l⋅s⎞             ⎛l⋅s⎞             ⎛l⋅s⎞     2     ⎛l⋅s⎞
   Rₗ⋅Rₛ⋅sinh⎜───⎟ + Rₗ⋅Z₀⋅cosh⎜───⎟ + Rₛ⋅Z₀⋅cosh⎜───⎟ + Z₀ ⋅sinh⎜───⎟
             ⎝ c ⎠             ⎝ c ⎠             ⎝ c ⎠           ⎝ c ⎠


The Laplace domain voltage across the input of the transmission line is

   >>> a.P1.V(s)
                           ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞⎞
                     Vₛ⋅Z₀⋅⎜Rₗ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟⎟
                           ⎝       ⎝ c ⎠          ⎝ c ⎠⎠
   ───────────────────────────────────────────────────────────────────────
     ⎛          ⎛l⋅s⎞             ⎛l⋅s⎞             ⎛l⋅s⎞     2     ⎛l⋅s⎞⎞
   s⋅⎜Rₗ⋅Rₛ⋅sinh⎜───⎟ + Rₗ⋅Z₀⋅cosh⎜───⎟ + Rₛ⋅Z₀⋅cosh⎜───⎟ + Z₀ ⋅sinh⎜───⎟⎟
     ⎝          ⎝ c ⎠             ⎝ c ⎠             ⎝ c ⎠           ⎝ c ⎠⎠

and the Laplace domain voltage across the output of the transmission line is

   >>> a.P2.V(s)
                                  Rₗ⋅Vₛ⋅Z₀
   ───────────────────────────────────────────────────────────────────────
     ⎛          ⎛l⋅s⎞             ⎛l⋅s⎞             ⎛l⋅s⎞     2     ⎛l⋅s⎞⎞
   s⋅⎜Rₗ⋅Rₛ⋅sinh⎜───⎟ + Rₗ⋅Z₀⋅cosh⎜───⎟ + Rₛ⋅Z₀⋅cosh⎜───⎟ + Z₀ ⋅sinh⎜───⎟⎟
     ⎝          ⎝ c ⎠             ⎝ c ⎠             ⎝ c ⎠           ⎝ c ⎠⎠


Since the transmission line is lossless, the transient response has an infinite number of reflections.  At the end of the transmission line

    >>> a.P2.V(t)
               ∞
             _____
             ╲
              ╲                                 m
               ╲   ⎛                          2⎞
                ╲  ⎜Rₗ⋅Rₛ - Rₗ⋅Z₀ - Rₛ⋅Z₀ + Z₀ ⎟   ⎛    l⋅(2⋅m + 1)⎞
    Rₗ⋅Vₛ⋅Z₀⋅   ╱  ⎜───────────────────────────⎟ ⋅u⎜t - ───────────⎟
               ╱   ⎜                          2⎟   ⎝         c     ⎠
              ╱    ⎝Rₗ⋅Rₛ + Rₗ⋅Z₀ + Rₛ⋅Z₀ + Z₀ ⎠
             ╱
             ‾‾‾‾‾
            m = 0
    ────────────────────────────────────────────────────────────────
                                                2
                      Rₗ⋅Rₛ + Rₗ⋅Z₀ + Rₛ⋅Z₀ + Z₀

The transient response at the start of the transmission line can be found in a similar way but is too long to show here.


If the transmission line is terminated at the load with its characteristic impedance:

   >>> b = a.subs({'Rl':'Z0'})
   >>> b.P1.V(s)
   ⎛ Vₛ⋅Z₀ ⎞
   ⎜───────⎟
   ⎝Rₛ + Z₀⎠
   ─────────
       s
   >>> b.P1.V(t)
   Vₛ⋅Z₀⋅u(t)
   ──────────
    Rₛ + Z₀
   >>> b.P2.V(s)
                                  Vₛ⋅Z₀
   ─────────────────────────────────────────────────────────────
     ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞⎞
   s⋅⎜Rₛ⋅sinh⎜───⎟ + Rₛ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟ + Z₀⋅cosh⎜───⎟⎟
     ⎝       ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠⎠
   >>> b.P2.V(t)
          ⎛    l⎞
   Vₛ⋅Z₀⋅u⎜t - ─⎟
          ⎝    c⎠
   ──────────────
      Rₛ + Z₀

If the transmission line is terminated at the source with its characteristic impedance:

   >>> c = a.subs({'Rs':'Z0'})
   >>> c.P1.V(s)
                     ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞⎞
                  Vₛ⋅⎜Rₗ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟⎟
                     ⎝       ⎝ c ⎠          ⎝ c ⎠⎠
   ─────────────────────────────────────────────────────────────
     ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞⎞
   s⋅⎜Rₗ⋅sinh⎜───⎟ + Rₗ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟ + Z₀⋅cosh⎜───⎟⎟
     ⎝       ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠⎠
   >>> c.P1.V(t)
      ⎛    2     2⎞ ⎛u(t)    ⎛    2⋅l⋅m⎞⎞
   Vₛ⋅⎝- Rₗ  + Z₀ ⎠⋅⎜──── + u⎜t - ─────⎟⎟
                    ⎝ 2      ⎝      c  ⎠⎠
   ──────────────────────────────────────
        -Rₗ⋅(Rₗ + Z₀) + Z₀⋅(Rₗ + Z₀)
   >>> c.P1.V(t).simplify()
      ⎛u(t)    ⎛    2⋅l⋅m⎞⎞
   Vₛ⋅⎜──── + u⎜t - ─────⎟⎟
      ⎝ 2      ⎝      c  ⎠⎠
   >>> c.P2.V(s)
                               Rₗ⋅Vₛ
   ─────────────────────────────────────────────────────────────
     ⎛       ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞          ⎛l⋅s⎞⎞
   s⋅⎜Rₗ⋅sinh⎜───⎟ + Rₗ⋅cosh⎜───⎟ + Z₀⋅sinh⎜───⎟ + Z₀⋅cosh⎜───⎟⎟
     ⎝       ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠          ⎝ c ⎠⎠
   >>> c.P2.V(t)
          ⎛    l⎞
   Rₗ⋅Vₛ⋅u⎜t - ─⎟
          ⎝    c⎠
   ──────────────
      Rₗ + Z₀

In this case, there is a single reflection.



Piezoelectric transducers
=========================

A piezoelectric transducer can be modelled with the KLM model.  This
can be described with the following netlist:

.. literalinclude:: examples/tutorials/transducers/KLM2.sch

Lcapy can draw this using:

   >>> from lcapy import Circuit, s, t
   >>> a = Circuit('KLM2.sch')
   >>> a.draw()

.. image:: examples/tutorials/transducers/KLM2.png
   :width: 12cm

In the Laplace domain, the transformer coupling ratio is

:math:`\phi = \frac{s Z_0}{2 h} \frac{1}{\sinh\left(\frac{s d}{2 c}\right)}`

and the impedance of `X1` is

:math:`X_1 = \frac{j h^2}{s^2 Z_0} \sinh\left(\frac{s d}{c}\right)`

where :math:`Z_0` is the characteristic impedance of the pizeoelectric
crystal, :math:`h` is the pressure constant of the crystal, and
:math:`d` is the thickness of the crysal.

The model has three ports: an electrical port, a back mechanical port,
and a fron mechanical port.


Transformers
============

The Lcapy `TF` component represents an ideal transformer.  Indeed it
operates at DC!


To model a real transformer, a pair of magnetically coupled inductors
is required.  The coupling is specified by the `K` component.  This
has three arguments: the names of two inductors and a coupling
coefficient `k`.   The mutual inductance of the two inductors is:

:math:`M = k \sqrt{L_1 L_2}`.

The inductances are related by the number of turns by:

:math:`L = \frac{N}{\mathcal{R}}`

where :math:`\mathcal{R}` is the reluctance.  If both inductors have a similar reluctance:

:math:`\frac{N_1}{N_2} = \sqrt{\frac{L_1}{L_2}}`

Note, this is an approximation and depends on the inductor geometry and the cores.

.. image:: examples/tutorials/transformers/mutualinductances1.png
   :width: 8cm

The equations for a pair of coupled inductors is

:math:`v_1(t) = L_1 \frac{\mathrm{d}}{\mathrm{d} t} i_1(t) + M \frac{\mathrm{d}}{\mathrm{d} t} i_2(t)`

:math:`v_2(t) = M \frac{\mathrm{d}}{\mathrm{d} t} i_1(t) + L_2 \frac{\mathrm{d}}{\mathrm{d} t} i_2(t)`

where :math:`M` is the mutual inductance.  This circuit can be
described in Lcapy with the netlist:

.. literalinclude:: examples/tutorials/transformers/mutualinductances1.sch

However, it cannot be analysed since the secondary has no reference to
the primary.  One approach is to couple the primary and secondary
circuits with a resistor and consider the result in the limit as the
resistor is made infinite.  There is no problem if the primary and
secondary share a common ground as in the figure below:

.. image:: examples/tutorials/transformers/mutualinductances2.png
   :width: 8cm


A transformer with a common ground for the primary and seconday can be
represented by a tee-model without having to use a coupling `K`
component.  In this case, the voltages can be expressed in terms of
the currents using:

:math:`v_1(t) = \left(L_1 - M\right) \frac{\mathrm{d}}{\mathrm{d} t} i_1(t) + M \frac{\mathrm{d}}{\mathrm{d} t} \left(i_1 - i_2(t)\right)`

:math:`v_2(t) = M \frac{\mathrm{d}}{\mathrm{d} t} \left(i_1(t) - i_2(t)\right) + \left(L_2 - M\right) \frac{\mathrm{d}}{\mathrm{d} t} i_2(t)`

Here :math:`L_1` is the primary magnetising inductance, :math:`L_2` is
the secondary magnetising inductance, :math:`L_1 - M` is the primary
leakage inductance, and :math:`L_2 - M` is the secondary leakage
inductance.  These equations correspond to the following circuit:

.. image:: examples/tutorials/transformers/tee-model-transformer1.png
   :width: 6.5cm


Sometimes a transformer is modelled around an ideal transformer with turns ratio :math:`a = N_1 / N_2` as per the following circuit:

.. image:: examples/tutorials/transformers/tee-model-transformer2.png
   :width: 12cm

Sometimes the secondary leakage inductance is referred to the primary:

.. image:: examples/tutorials/transformers/tee-model-transformer3.png
   :width: 11cm
