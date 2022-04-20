S1 box; right=1.9, aspect=2.1, draw=white, l={Laplace domain\\$H(s)$}
S2 box; right=1.9, aspect=2.1, draw=white, l={Time domain\\$h(t)$}
S3 box; right=1.9, aspect=2.1, draw=white, l={Phasor domain\\$H(\mathrm{j}\omega)$}
S4 box; right, l={$\mathcal{L}\left\{.\right\}$}
S5 box; right=1.9, aspect=2.1, draw=white, l={Angular Fourier domain\\$H(\omega)$}
S6 box; right, l={$s=\mathrm{j}\omega$}
S7 box; right, l={$\mathcal{F}_{\omega}\left\{.\right\}$}
S8 box; right=1.9, aspect=2.1, draw=white, l={Equivalent if lossy and causal}
S9 box; right, l={$\mathcal{F}\left\{.\right\}$}
S10 box; right=1.9, aspect=2.1, draw=white, l={Fourier domain\\$H(f)$}
S11 box; right, l={$f=\frac{\omega}{2\pi}$}

W S1.e S4.w; right=0.4, startarrow=tri, color=blue
W S4.e S2.w; right=0.4, startarrow=tri, color=black!70!green
W S2.s S7.n; down=0.4, endarrow=tri, color=black!70!green
W S7.s S5.n; down=0.4, endarrow=tri, color=green
W S1.s S6.n; down=0.4, endarrow=tri, color=blue
W S6.s S3.n; down=0.4, endarrow=tri, color=purple
W S3.e S8.w; right=0.4, dashed, color=purple
W S8.e S5.w; right=0.4, dashed, color=green

W S2.e S9.w; right=0.4, endarrow=tri, color=black!70!green
W S9.e S10.w; right=0.4, endarrow=tri, color=magenta
W S10.s S11.n; down=0.4, endarrow=tri, color=magenta
W S11.s 5; down=0.4, color=green
W 5 S5.e; left, endarrow=tri, color=green

; label_nodes=alpha, draw_nodes=none, label_ids=false
