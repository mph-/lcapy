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
\usetikzlibrary{fit, shapes, arrows, patterns, decorations.text, decorations.markings}

\begin{document}
\begin{tikzpicture}[scale=1.00, transform shape, /tikz/circuitikz/bipoles/length=1.50cm, american currents, american voltages, voltage dir=RP]
  \coordinate (1) at (0,0);
  \coordinate (2) at (2,0);
  \coordinate (3) at (2.5,0);
  \coordinate (0) at (4.5,0);
  \draw (1) to [R, l_=$1\,\mbox{$\Omega$}$, n=R1] (2);
  \draw[-] (2) to (3);
  \draw (3) to [C, l_=$2\,\mbox{F}$, n=C1] (0);
  \draw (1) node[ocirc] {};
  \draw (0) node[ocirc] {};
\end{tikzpicture}
\end{document}"""

        content = open('tmp.tex').read()
        if content != match:
            raise ValueError('Schematic mismatch')

    def test_schematic2(self):

        a = Circuit("""
R1 1 2 1; right
R2 2 3 2; right
C 1 3 3; right""")
        a.draw('tmp.tex')

        match = r"""\documentclass[a4paper]{standalone}
\usepackage{amsmath}
\usepackage{circuitikz}
\usetikzlibrary{fit, shapes, arrows, patterns, decorations.text, decorations.markings}

\begin{document}
\begin{tikzpicture}[scale=1.00, transform shape, /tikz/circuitikz/bipoles/length=1.50cm, american currents, american voltages, voltage dir=RP]
  \coordinate (1) at (0,0);
  \coordinate (2) at (2,0);
  \coordinate (3) at (4,0);
  \draw (1) to [R, l_={$R_{1}$}{=$1\,\mbox{$\Omega$}$}, n=R1] (2);
  \draw (2) to [R, l_={$R_{2}$}{=$2\,\mbox{$\Omega$}$}, n=R2] (3);
  \draw (1) to [C, l_={$C$}{=$3\,\mbox{F}$}, n=C] (3);
  \draw (1) node[circ] {};
  \draw (2) node[circ] {};
  \draw (3) node[circ] {};
  \draw[anchor=south east] (1) node {1};
  \draw[anchor=south east] (2) node {2};
  \draw[anchor=south east] (3) node {3};
\end{tikzpicture}
\end{document}"""

        content = open('tmp.tex').read()
        if content != match:
            raise ValueError('Schematic mismatch')

    def test_network_node_positions1(self):

        a = R(1) + C(2)
        sch = a.sch()
        sch._positions_calculate()
        nodes = sch.nodes

        self.assertEqual(nodes['1'].pos.x, 0, '1 x')
        self.assertEqual(nodes['2'].pos.x, 2, '2 x')
        self.assertEqual(nodes['3'].pos.x, 2.5, '3 x')
        self.assertEqual(nodes['0'].pos.x, 4.5, '0 x')

    def test_network_node_positions2(self):

        a = (R(1) + C(2)) | L(3)
        sch = a.sch()
        sch._positions_calculate()
        nodes = sch.nodes

        self.assertEqual(nodes['1'].pos.x, 0, '1 x')
        self.assertEqual(nodes['2'].pos.x, 0.5, '2 x')
        self.assertEqual(nodes['3'].pos.x, 5, '3 x')
        self.assertEqual(nodes['4'].pos.x, 0.5, '4 x')
        self.assertEqual(nodes['5'].pos.x, 5, '5 x')
        self.assertEqual(nodes['6'].pos.x, 2.5, '6 x')
        self.assertEqual(nodes['7'].pos.x, 3, '7 x')
        self.assertEqual(nodes['8'].pos.x, 0.5, '8 x')
        self.assertEqual(nodes['9'].pos.x, 5, '9 x')
        self.assertEqual(nodes['0'].pos.x, 5.5, '0 x')

    def test_network_danglin3(self):

        a = Circuit("""
R1 1 2; right
R2 2 3; right
R3 2 5; down
R4 4 5; right=2
R5 5 6; right
R6 5 7; down
R7 7 8; right=2, fixed""")

        sch = a.sch
        sch._positions_calculate()
        nodes = sch.nodes

        self.assertEqual(nodes['1'].pos.x, 2, '1 x')
        self.assertEqual(nodes['2'].pos.x, 4, '2 x')
        self.assertEqual(nodes['3'].pos.x, 6, '3 x')
        self.assertEqual(nodes['4'].pos.x, 0, '4 x')
        self.assertEqual(nodes['5'].pos.x, 4, '5 x')
        self.assertEqual(nodes['6'].pos.x, 6, '6 x')
        self.assertEqual(nodes['7'].pos.x, 4, '7 x')
        self.assertEqual(nodes['8'].pos.x, 8, '8 x')
