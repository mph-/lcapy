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


print("""
.. _schematic gallery:

==================
Schematic gallery
==================
""")


# TODO: add some curation...
filenames = glob('examples/**/*.sch', recursive=True)
filenames.sort(key=lambda x: basename(x))

for filename in filenames:
    show_image(filename)
