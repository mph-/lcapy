=====================
Discrete-time signals
=====================


The discrete-time signal support is experimental and probably riddled
with bugs.  Importing the `discretime` module, introduces three new
domain variables:

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

The z-transform is performed explicitly with the `ZT` method:

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


Z-transform expressions are objects of the `zExpr` class.  They are functions of the complex variable `z` and are similar to `sExpr` objects.

The poles and zeros can be plotted using the `plot()` method.  For example,

.. literalinclude:: examples/discretetime/dt1-pole-zero-plot1.py

.. image:: examples/discretetime/dt1-pole-zero-plot1.png
   :width: 15cm


Discrete time Fourier transform (DTFT)
======================================

The DTFT converts an n-domain or z-domain expression into the f-domain (continuous Fourier domain).  For example,

.. literalinclude:: examples/discretetime/dt1-DTFT-plot1.py

.. image:: examples/discretetime/dt1-DTFT-plot1.png
   :width: 15cm


