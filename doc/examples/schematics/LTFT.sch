S1 box; right=1.9, aspect=2.1, draw=white, l={Transfer function\\$H(s)$}
S2 box; right=1.9, aspect=2.1, draw=white, l={Impulse response\\$h(t)$}
S3 box; right=1.9, aspect=2.1, draw=white, l={Frequency response\\$H(\mathrm{j}2\pi f)$ or $H(f)$}
S4 box; right, l={$\mathcal{L}^{-1}\left\{.\right\}$}
S5 box; right, l={$\mathcal{F}\left\{.\right\}$}
S6 box; right, l={$s=\mathrm{j}2\pi f$}
W S1.e S4.w; right=0.4, startarrow=tri, endarrow=tri, color=blue
W S4.e S2.w; right=0.4, startarrow=tri, endarrow=tri, color=black!70!green
W S2.s 4; down, color=black!70!green
W 4 S5.e; left=0.4, endarrow=tri, color=black!70!green
W S5.w S3.e; left=0.4, startarrow=tri, endarrow=tri, color=purple
W S1.s S6.n; down=0.4, endarrow=tri, color=blue
W S6.s S3.n; down=0.4, endarrow=tri, color=purple
O S4.mid S5.mid; down=1.3
; label_nodes=alpha, draw_nodes=none, label_ids=false

