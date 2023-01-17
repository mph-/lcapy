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
    'E1', 'F1', 'G1', 'H1', 'K2', 'VRL2',
    'VRmesh1',
    'cmos-esd-damage', 'cmos-input-model1', 'cmos-led1', 'cmos-open-drain',
    'cmos-protection1', 'cmos-totem2', 'cmos_C_load_thevenin',
    'cmos_R_series_C_load_thevenin', 'common-base',
    'ferrite-choke1', 'filter1', 'filter2', 'fit3',
    'non-ideal-transformer-primary', 'non-ideal-transformer')))

categories.append(Category('Block diagrams', ('LTFT', 'domains',
                  'fir5', 'negative-feedback2', 'transforms')))

categories.append(Category('Two port networks', (
    'TPA', 'TPB', 'TPG', 'TPH', 'TPY', 'TPZ', 'series', 'shunt',
    'pisection', 'tsection')))

categories.append(Category('Components', (
    'ADC1', 'arrows', 'antennas', 'batteries', 'buffers',
    'cable-coax', 'cable-tp', 'cable-twinax', 'cable-utp', 'capacitors',
    'connections1', 'connections2', 'CPE1', 'DAC3', 'diodes', 'diodes2',
    'FB1', 'fdopamps', 'flipflops', 'grounds', 'GY1',
    'hsources', 'massspringdamper1', 'meters1',
    'misc', 'multiplexers', 'opamps',
    'resistors1', 'resistors2', 'resistors3', 'resistors4',
    'resistors5', 'resistors6', 'resistors7', 'switches', 'switches2',
    'transformers', 'transistors', 'transistors2',
    'Uinamp', 'Uisoamp', 'Ufdopamp', 'Uopamp',
    'variable1', 'vsources', 'wirestyles', 'XT1')))

categories.append(Category('Chips', ('chips', 'chips2', 'chips3')))

categories.append(
    Category('Shapes', ('Striangle2', 'Sbox1', 'Scircle2', 'circles')))

categories.append(Category('Labels', (
    'annotate1', 'current_labels1', 'current_labels2',
    'flow_labels1', 'labels1',
    'voltage_labels1', 'Rlabels')))

categories.append(Category('Digital circuits', (
    'PIO1', 'PIO-input2', 'PIO-input3', 'PIO-output2',  'PIO-pullup2', 'mcu-bt',
    'cmos-backdrive1', 'cmos-backdrive2', 'I2C2', 'ic1',
    'stepup')))

categories.append(Category('Transistor circuits',
                  ('nmosfet_simple_amplifier', 'mosfet_push_pull_load',
                   'BH-mag')))

categories.append(Category('Opamp circuits', (
    'opamp-differential-amplifier1',
    'opamp-inverting-amplifier1',
    'opamp-inverting-integrator',
    'opamp-transimpedance-amplifier1',
    'opamp-transimpedance-amplifier2',
    'inamp-amplifier1', 'fdopamp-amplifier1',
    'opamp-displacement-current-sensor-noise-model1',
    'opamp-displacement-current-sensor1',
    'opamp-noninverting-amplifier-noisy', 'opamp-piezo-amplifier1',
    'opamp-transimpedance-amplifier-with-voltage-gain1',
    'pierce-oscillator', 'sallen-key-lpf1', 'multiple-feedback-lpf',
    'guard1', 'shield-ground', 'shield-guard',
    'opamp-displacement-current-sensor-noise-model1',
    'opamp-displacement-current-sensor1',
    'opamp-noninverting-amplifier-noisy')))


categories.append(Category('Transducers',
                           ('KLM2', 'KLM3', 'accelerometer1')))

categories.append(Category('Miscellaneous', ('colors1', )))

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

The images in this section are bit-mapped PNG images with DPI=300.  Higher quality images can be generated using PDF or PGF formats.

""")

for category in categories:
    category.show(filename_dict)
