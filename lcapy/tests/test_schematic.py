from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_schematic1(self):

        a = R(1) + C(2)
        a.draw('tmp.tex')

        match = r"""\documentclass[a4paper]{standalone}
\usepackage{amsmath}
\usepackage{circuitikz}
\usetikzlibrary{fit, shapes, arrows}
\begin{document}
\begin{tikzpicture}[scale=1.00, transform shape, /tikz/circuitikz/bipoles/length=1.50cm, american currents, american voltages]
  \coordinate (1) at (0,0);
  \coordinate (2) at (2,0);
  \coordinate (3) at (3,0);
  \coordinate (0) at (5,0);
  \draw[] (1) to [R,l_=$1\,\mbox{$\Omega$}$,,,,o-,n=R1] (2);
  \draw (1) node[ocirc] {};
  \draw[-, , ] (2) to (3);
  \draw[] (3) to [C,l_=$2\,\mbox{F}$,,,,-o,n=C1] (0);
  \draw (0) node[ocirc] {};
\end{tikzpicture}
\end{document}"""
        
        content = open('tmp.tex').read()
        if content != match:
            raise Error('Schematic mismatch')


    def test_schematic2(self):

        a = Circuit("""
R1 1 2 1; right
R2 2 3 2; right
C 1 3 3; right""")
        a.draw('tmp.tex')

        match = r"""\documentclass[a4paper]{standalone}
\usepackage{amsmath}
\usepackage{circuitikz}
\usetikzlibrary{fit, shapes, arrows}
\begin{document}
\begin{tikzpicture}[scale=1.00, transform shape, /tikz/circuitikz/bipoles/length=1.50cm, american currents, american voltages]
  \coordinate (1) at (0,0);
  \coordinate (2) at (2,0);
  \coordinate (3) at (4,0);
  \draw[] (1) to [R,l_={$R_{1}$}{=$1\,\mbox{$\Omega$}$},,,,*-*,n=R1] (2);
  \draw (1) node[circ] {};
  \draw (2) node[circ] {};
  \draw[] (2) to [R,l_={$R_{2}$}{=$2\,\mbox{$\Omega$}$},,,,*-*,n=R2] (3);
  \draw (2) node[circ] {};
  \draw (3) node[circ] {};
  \draw[] (1) to [C,l_={C}{=$3\,\mbox{F}$},,,,*-*,n=C] (3);
  \draw (1) node[circ] {};
  \draw (3) node[circ] {};
  \draw[anchor=south east] (1) node {1};
  \draw[anchor=south east] (2) node {2};
  \draw[anchor=south east] (2) node {2};
  \draw[anchor=south east] (3) node {3};
  \draw[anchor=south east] (1) node {1};
  \draw[anchor=south east] (3) node {3};
\end{tikzpicture}
\end{document}"""
        
        content = open('tmp.tex').read()
        if content != match:
            raise Error('Schematic mismatch')
        
