=====================
Discrete-time signals
=====================

Discrete-time signal support is experimental and requires the  `discretetime` to be explicitly imported.  It introduces three new domain variables:

  - `n` for discrete-time signals, for example, `3 * u(n - 2)`
  - `k` for discrete-frequency spectra
  - `z` for z-transforms, for example, `Y(z)`

The `n`, `k`, and `z` variables share many of the attributes and methods of their continuous-time equivalents, `t`, `f`, and `s`, :ref:`expressions`.
    

The discrete-time signal can be plotted using the `plot()` method.
For example,

.. literalinclude:: examples/discretetime/dt1-plot1.py

.. image:: examples/discretetime/dt1-plot1.png
   :width: 15cm   
   

Sequences
=========

   >>> x = unitimpulse(n) + 2 * unitimpulse(n - 2)
   >>> seq = x.seq((-5, 5))
   >>> seq
       {0, 0, 0, 0, 0, _1, 0, 2, 0, 0, 0}

Note, the underscore marks the item in the sequence where `n = 0`.
       
   >>> seq.as_impulses()
   >>> δ[n] + 2⋅δ[n - 2]

   >>> seq.extent()
   >>> 3


Z-transform
===========

Lcapy uses the unilateral z-transform.  It is performed explicitly with the `ZT` method:

   >>> x = unitimpulse(n) + 2 * unitimpulse(n - 2)
   >>> x.ZT()
   >>>      2 
       1 + ──
            2
           z

It is also performed implicitly with `z` as an argument:
      
   >>> x(z)
   >>>     2 
      1 + ──
           2
          z

The inverse unilateral z-transform is not unique and is only defined for :math:`n \ge 0`.  For example,

   >>> H = z / (z - 'a')
   >>> H(n)
   ⎧ n           
   ⎨a   for n ≥ 0
   ⎩             

   
If the result is known to be causal, then use:   
   >>> H(n, causal=True)
    n     
   a ⋅u(n)


Z-domain expressions are objects of the `zExpr` class.  They are functions of the complex variable `z` and are similar to `sExpr` objects.


The poles and zeros of a z-domain expression can be plotted using the `plot()` method.  For example,

.. literalinclude:: examples/discretetime/dt1-pole-zero-plot1.py

.. image:: examples/discretetime/dt1-pole-zero-plot1.png
   :width: 15cm


Discrete time Fourier transform (DTFT)
======================================

The DTFT converts an n-domain or z-domain expression into the f-domain (continuous Fourier domain).  Note, unlike the Fourier transform, this is periodic with period :math:`1/\Delta t`.  For example,

.. literalinclude:: examples/discretetime/dt1-DTFT-plot1.py

.. image:: examples/discretetime/dt1-DTFT-plot1.png
   :width: 15cm


Bilinear transform
==================

The bilinear transform can be used to approximate an s-domain expression with a z-domain expression using :math:`s \approx \frac{2}{\Delta t} \frac{1 - z^{-1}}{1 + z^{-1}}`.   This is performed by the `bilinear_transform()` method of s-domain objects, for example,

   >>> H = s / (s - 'a')
   >>> Hz = H.bilinear_transform().simplify()
   >>> Hz
         2⋅(1 - z)       
   ──────────────────────
   Δₜ⋅a⋅(z + 1) - 2⋅z + 2

