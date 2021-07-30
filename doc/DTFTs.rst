:math:`x(n) \longleftrightarrow \sum_{m=-\infty}^{\infty} X(f - \frac{m}{\Delta_{t}})`

:math:`x(a n) \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} X(\frac{f}{a} - \frac{m}{\Delta_{t} a})}{\left|{a}\right|}`

:math:`x(- m + n) \longleftrightarrow \sum_{p=-\infty}^{\infty} X(f - \frac{p}{\Delta_{t}}) e^{- 2 \mathrm{j} \pi f m} e^{\frac{2 \mathrm{j} \pi m p}{\Delta_{t}}}`

:math:`\cos{\left(2 \pi \Delta_{t} f_{0} n \right)} \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} \delta\left(f - f_{0} - \frac{m}{\Delta_{t}}\right)}{2 \Delta_{t}} + \frac{\sum_{m=-\infty}^{\infty} \delta\left(f + f_{0} - \frac{m}{\Delta_{t}}\right)}{2 \Delta_{t}}`

:math:`\sin{\left(2 \pi \Delta_{t} f_{0} n \right)} \longleftrightarrow - \frac{\mathrm{j} \sum_{m=-\infty}^{\infty} \delta\left(f - f_{0} - \frac{m}{\Delta_{t}}\right)}{2 \Delta_{t}} + \frac{\mathrm{j} \sum_{m=-\infty}^{\infty} \delta\left(f + f_{0} - \frac{m}{\Delta_{t}}\right)}{2 \Delta_{t}}`

:math:`e^{2 \mathrm{j} \pi \Delta_{t} f_{0} n} \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} \delta\left(f - f_{0} - \frac{m}{\Delta_{t}}\right)}{\Delta_{t}}`

:math:`1 \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} \delta\left(f - \frac{m}{\Delta_{t}}\right)}{\Delta_{t}}`

:math:`\delta\left[n\right] \longleftrightarrow 1`

:math:`\delta\left[- m + n\right] \longleftrightarrow e^{- 2 \mathrm{j} \pi \Delta_{t} f m}`

:math:`u\left[n\right] \longleftrightarrow \frac{e^{2 \mathrm{j} \pi \Delta_{t} f}}{e^{2 \mathrm{j} \pi \Delta_{t} f} - 1} + \frac{\sum_{m=-\infty}^{\infty} \delta\left(f - \frac{m}{\Delta_{t}}\right)}{2 \Delta_{t}}`

:math:`n u\left[n\right] \longleftrightarrow \frac{e^{2 \mathrm{j} \pi \Delta_{t} f}}{e^{4 \mathrm{j} \pi \Delta_{t} f} - 2 e^{2 \mathrm{j} \pi \Delta_{t} f} + 1} + \frac{\mathrm{j} \sum_{m=-\infty}^{\infty} \delta^{\left( 1 \right)}\left( f - \frac{m}{\Delta_{t}} \right)}{4 \pi \Delta_{t}^{2}}`

:math:`\mathop{\mathrm{sign}}\left[n\right] \longleftrightarrow \frac{2 e^{2 \mathrm{j} \pi \Delta_{t} f}}{e^{2 \mathrm{j} \pi \Delta_{t} f} - 1}`

:math:`\alpha^{- n} u\left[n\right] \longleftrightarrow \frac{\alpha e^{2 \mathrm{j} \pi \Delta_{t} f}}{\alpha e^{2 \mathrm{j} \pi \Delta_{t} f} - 1}`

:math:`\mathop{\mathrm{rect}}\left[n\right] \longleftrightarrow 1`

:math:`\mathop{\mathrm{rect}}\left[\frac{n}{N_{o}}\right] \longleftrightarrow \frac{\sin{\left(\pi \Delta_{t} N_{o} f \right)}}{\sin{\left(\pi \Delta_{t} f \right)}}`

:math:`\mathop{\mathrm{rect}}\left[\frac{n}{N_{e}}\right] \longleftrightarrow \frac{e^{\mathrm{j} \pi \Delta_{t} f} \sin{\left(\pi \Delta_{t} N_{e} f \right)}}{\sin{\left(\pi \Delta_{t} f \right)}}`

:math:`\mathrm{sincn}{\left(n \right)} \longleftrightarrow \sum_{m=-\infty}^{\infty} \mathop{\mathrm{rect}}\left[\Delta_{t} f - m\right]`

:math:`\mathrm{sincn}{\left(K n \right)} \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} \mathop{\mathrm{rect}}\left[\frac{\Delta_{t} f}{K} - \frac{m}{K}\right]}{K}`

:math:`\mathrm{sincn}^{2}{\left(K n \right)} \longleftrightarrow \frac{\sum_{m=-\infty}^{\infty} \operatorname{tri}{\left(\frac{\Delta_{t} f}{K} - \frac{m}{K} \right)}}{K}`

:math:`\mathrm{sincu}{\left(K n \right)} \longleftrightarrow \frac{\pi \sum_{m=-\infty}^{\infty} \mathop{\mathrm{rect}}\left[\frac{\pi \Delta_{t} f}{K} - \frac{\pi m}{K}\right]}{K}`

