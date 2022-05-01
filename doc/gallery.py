import pdb
from struct import unpack
from glob import glob
from os.path import basename
import imghdr


def get_image_size(filename):

    with open(filename, 'rb') as fhandle:
        header = fhandle.read(24)
        if len(header) != 24:
            raise ValueError('Truncated file')
        if imghdr.what(filename) == 'png':
            check = unpack('>i', header[4:8])[0]
            if check != 0x0d0a1a0a:
                raise ValueError('Not PNG file')
            width, height = unpack('>ii', header[16:24])
    return width, height


def show_image(sch_filename):

    png_filename = sch_filename.replace('.sch', '.png')
    try:
        width, height = get_image_size(png_filename)
    except:
        return

    width_cm = width / 120

    label = basename(sch_filename)

    url = r'https://github.com/mph-/lcapy/tree/master/doc/' + sch_filename

    link = "`%s <%s>`_" % (label, url)

    print(
        """
%s

.. literalinclude:: %s

.. image:: %s
   :width: %.1fcm
        """ % (link, sch_filename, png_filename, width_cm))


class Category:

    def __init__(self, name, basenames):

        self.name = name
        self.basenames = basenames

    def show(self, filename_dict):

        header = self.name
        underline = '=' * len(self.name)
        print(header)
        print(underline + '\n')

        for name in self.basenames:
            filename = filename_dict[name + '.sch']
            show_image(filename)


categories = []

categories.append(Category('Circuits', (
    'E1', 'F1', 'G1', 'H1', 'I2C2', 'K2', 'KLM2', 'KLM3', 'VRL2',
    'VRmesh1', 'accelerometer1', 'BH-mag', 'cmos-backdrive1', 'cmos-backdrive2',
    'cmos-esd-damage', 'cmos-input-model1', 'cmos-led1', 'cmos-open-drain',
    'cmos-protection1', 'cmos-totem2', 'cmos_C_load_thevenin',
    'cmos_R_series_C_load_thevenin', 'common-base', 'fdopamp-amplifier1',
    'ferrite-choke1', 'filter1', 'filter2', 'fit3', 'guard1', 'ic1',
    'inamp-amplifier1', 'integrator', 'mcu-bt', 'multiple-feedback-lpf',
    'non-ideal-transformer-primary', 'non-ideal-transformer',
    'opamp-differential-amplifier1',
    'opamp-displacement-current-sensor-noise-model1',
    'opamp-displacement-current-sensor1', 'opamp-inverting-amplifier1',
    'opamp-inverting-integrator', 'opamp-noninverting-amplifier-noisy',
    'opamp-piezo-amplifier1',
    'opamp-transimpedance-amplifier-with-voltage-gain1',
    'opamp-transimpedance-amplifier1',
    'opamp-transimpedance-amplifier2', 'pierce-oscillator',
    'sallen-key-lpf1', 'shield-ground', 'shield-guard', 'stepup')))

categories.append(Category('Block diagrams', ('LTFT', 'domains',
                  'fir5', 'negative-feedback2', 'transforms')))

categories.append(Category('Two port networks', ('TPA', 'TPB', 'TPG',
                  'TPH', 'TPY', 'TPZ', 'series', 'shunt')))

categories.append(Category('Components', (
    'cable-tp', 'cable-twinax', 'cable-utp', 'connections1', 'diodes',
    'diodes2', 'fdopamps', 'grounds', 'hsources', 'massspringdamper1',
    'misc', 'opamps', 'resistors1', 'resistors2', 'resistors3', 'resistors4',
    'resistors5', 'resistors6', 'resistors7', 'switches', 'switches2',
    'transistors', 'transistors2', 'variable1', 'vsources', 'wirestyles')))

categories.append(Category('Chips', ('chips', 'chips2', 'chips3')))

categories.append(Category('Labels', ('current_labels1',
                  'current_labels2', 'flow_labels1', 'voltage_labels1')))

categories.append(Category('Miscellaneous', ('circles', 'colors1')))

filename_dict = {}
filenames = glob('examples/**/*.sch', recursive=True)
filenames.sort(key=lambda x: basename(x))
for filename in filenames:
    filename_dict[basename(filename)] = filename


print("""
.. _schematic gallery:

==================
Schematic gallery
==================
""")

for category in categories:
    category.show(filename_dict)
