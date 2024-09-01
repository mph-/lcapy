from .chip import Chip


class Uadc(Chip):
    """ADC"""

    pins = {'in': ('l', -0.5, 0),
            'in+': ('l', -0.4375, 0.125),
            'in-': ('l', -0.4375, -0.125),
            'vref-': ('l', -0.375, -0.25),
            'vref+': ('l', -0.375, 0.25),
            'avss': ('b', -0.1, -0.5),
            'vss': ('b', 0.1, -0.5),
            'dvss': ('b', 0.3, -0.5),
            'clk': ('r', 0.5, -0.25),
            'data': ('r', 0.5, 0),
            'fs': ('r', 0.5, 0.25),
            'dvdd': ('t', 0.3, 0.5),
            'vdd': ('t', 0.1, 0.5),
            'avdd': ('t', -0.1, 0.5)}

    pinlabels = {'vref-': 'VREF-', 'vref+': 'VREF+',
                 'vss': 'VSS', 'vdd': 'VDD',
                 'dvss': 'DVSS', 'dvdd': 'DVDD',
                 'avss': 'AVSS', 'avdd': 'AVDD',
                 'clk': '>', 'data': 'DATA', 'fs': 'FS'}

    @property
    def path(self):
        return ((-0.5, 0.0), (-0.25, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.25, 0.5))
