from .chip import Chip


class Udac(Chip):
    """DAC"""

    pins = {'out': ('r', 0.5, 0),
            'out+': ('r', 0.4375, 0.125),
            'out-': ('r', 0.4375, -0.125),
            'vref-': ('r', 0.375, -0.25),
            'vref+': ('r', 0.375, 0.25),
            'avss': ('b', 0.1, -0.5),
            'vss': ('b', -0.1, -0.5),
            'dvss': ('b', -0.3, -0.5),
            'clk': ('l', -0.5, -0.25),
            'data': ('l', -0.5, 0),
            'fs': ('l', -0.5, 0.25),
            'dvdd': ('t', -0.3, 0.5),
            'vdd': ('t', -0.1, 0.5),
            'avdd': ('t', 0.1, 0.5)}

    pinlabels = {'vref-': 'VREF-', 'vref+': 'VREF+',
                 'vss': 'VSS', 'vdd': 'VDD',
                 'dvss': 'DVSS', 'dvdd': 'DVDD',
                 'avss': 'AVSS', 'avdd': 'AVDD',
                 'clk': '>', 'data': 'DATA', 'fs': 'FS'}

    @property
    def path(self):
        return ((-0.5, -0.5), (0.25, -0.5), (0.5, 0), (0.25, 0.5), (-0.5, 0.5))
