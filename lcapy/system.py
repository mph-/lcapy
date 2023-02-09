"""This module wraps system dependent programs for pdf generation, etc.

Copyright 2014--2021 Michael Hayes, UCECE

"""


from os import path, remove, chdir, getcwd, stat, system
from subprocess import call
import re
import platform

try:
    from subprocess import DEVNULL
except ImportError:
    # For python 2.7
    import os
    DEVNULL = open(os.devnull, 'wb')


# System dependent functions

# Note, sometime in 2018 ImageMagick's convert program by default
# disallowed pdf file conversions.  A work-around is to edit
# /etc/ImageMagick-6/policy.xml to enable this conversion.

programs = {}


def which(program):

    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    for path1 in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path1, program)
        if is_exe(exe_file):
            return exe_file

    return None


def hasexe(program):

    try:
        return programs[program]
    except:
        if which(program) is not None:
            programs[program] = True
        else:
            programs[program] = which(program + '.exe') is not None

    return programs[program]


def checkexe(command, debug=False):

    if isinstance(command, list):
        program = command[0]
        command = ' '.join(command)
    else:
        program = command.split(' ')[0]

    if hasexe(program):
        if debug:
            print('Checking: %s found' % program)
    else:
        if debug:
            print('Checking: %s not found' % program)
        raise RuntimeError('%s is not installed' % program)


def run(command, stderr=DEVNULL, stdout=DEVNULL, shell=False, debug=False):

    checkexe(command, debug=debug)
    if debug:
        print('Running: %s' % ' '.join(command))
    call(command, stderr=stderr, stdout=stdout, shell=shell)


def tmpfilename(suffix='', dirname=None):

    from tempfile import gettempdir, NamedTemporaryFile

    if dirname is None:
        # Searches using TMPDIR, TEMP, TMP environment variables
        dirname = gettempdir()

    filename = NamedTemporaryFile(suffix=suffix, dir=dirname,
                                  delete=False).name
    return filename


class PDFConverter(object):

    def __init__(self, debug=False):

        self.debug = debug

    def to_svg(self, pdf_filename, svg_filename):

        run(['pdf2svg', pdf_filename, svg_filename], debug=self.debug)
        if not path.exists(svg_filename):
            raise RuntimeError('Could not generate %s with pdf2svg.  Is it installed?' %
                               svg_filename)

    def to_png_convert(self, pdf_filename, png_filename, dpi=300):

        program = 'convert'
        if platform.system() == 'Windows':
            program = 'magick convert'

        run([program, '-density %d' %
            int(dpi), pdf_filename, png_filename], debug=self.debug)

        if not path.exists(png_filename):
            raise RuntimeError('Could not generate %s with convert' %
                               png_filename)
        if stat(png_filename).st_size == 0:
            raise RuntimeError('Could not generate %s with convert, empty file' %
                               png_filename)

    def to_png_ghostscript(self, pdf_filename, png_filename, dpi=300):

        program = 'gs'
        if platform.system() == 'Windows':
            program = 'gswin32'
            if platform.machine().endswith('64'):
                program = 'gswin64'

        run([program, '-q', '-dQUIET', '-dSAFER', '-dBATCH', '-dNOPAUSE',
             '-dNOPROMPT',  '-dMaxBitmap=500000000',  '-dAlignToPixels=0',
             '-dGridFitTT=2',  '-sDEVICE=pngalpha',  '-dTextAlphaBits=4',
             '-dGraphicsAlphaBits=4',  '-r%dx%d' % (int(dpi), int(dpi)),
             '-sOutputFile=' + png_filename, pdf_filename], debug=self.debug)

    def to_png_pdftoppm(self, pdf_filename, png_filename, dpi=300):

        root, ext = path.splitext(png_filename)
        args = ['pdftoppm', '-r %d' %
                int(dpi), '-png', '-thinlinemode shape', '-singlefile', pdf_filename, root]
        if False:
            # TODO, determine why this fails...
            run(args, debug=self.debug)
        else:
            checkexe(args)
            system(' '.join(args))

        if not path.exists(png_filename):
            raise RuntimeError('Could not generate %s with pdftoppm' %
                               png_filename)

    def to_png(self, pdf_filename, png_filename, dpi=300, method=None):

        converters = {'ghostscript': self.to_png_ghostscript,
                      'convert': self.to_png_convert,
                      'pdftoppm': self.to_png_pdftoppm}

        # Ordered in preference
        methods = ('ghostscript', 'convert', 'pdftoppm')

        if method is not None:
            if method not in methods:
                raise KeyError('Unknown method %s, known methods: %s' %
                               (method, ', '.join(methods)))

            # Try the specified method then fall back to default methods
            # if it does not work
            func = converters[method]
            return func(pdf_filename, png_filename, dpi)

        for method in methods:
            func = converters[method]
            try:
                return func(pdf_filename, png_filename, dpi)
            except:
                pass

        raise RuntimeError(
            """Could not convert pdf to png, tried: %s.  Check that one of these programs is installed.  See https://lcapy.readthedocs.io/en/latest/install.html""" % ', '.join(methods))


def run_dot(dotfilename, filename):

    base, ext = path.splitext(filename)
    #run([which('dot'), '-T ' + ext[1:], '-o ' + filename, dotfilename])
    system('dot -T %s -o %s %s' % (ext[1:], filename, dotfilename))
    remove(dotfilename)
    return


class LatexRunner(object):

    def __init__(self, debug=False):

        self.debug = debug

    def run(self, tex_filename):

        checkexe('pdflatex', debug=self.debug)

        root, ext = path.splitext(tex_filename)
        dirname = path.dirname(tex_filename)
        baseroot = path.basename(root)
        cwd = getcwd()
        if dirname != '':
            if self.debug:
                print('Chdir: %s' % dirname)
            chdir(path.abspath(dirname))

        run(['pdflatex', '-interaction', 'batchmode',
            baseroot + '.tex'], debug=self.debug)

        if dirname != '':
            if self.debug:
                print('Chdir: %s' % cwd)
            chdir(cwd)

        return root + '.pdf'

    def cleanup(self, tex_filename, wanted_filename=''):

        root, ext = path.splitext(tex_filename)
        exts = ['.tex', '.aux', '.log', '.pdf']
        filenames = [root + ext for ext in exts]

        if wanted_filename in filenames:
            filenames.remove(wanted_filename)

        for filename in filenames:
            if path.exists(filename):
                if self.debug:
                    print('Removing: %s' % filename)
                remove(filename)

    def extract_circuitikz_version(self, tex_filename):

        log_filename = tex_filename.replace('.tex', '.log')

        with open(log_filename, 'rt') as logfile:
            lines = logfile.readlines()

        date = None
        version = None
        for m, line in enumerate(lines):
            if line.startswith('Package: circuitikzgit'):
                match = re.search(r'circuitikzgit ([0-9/]+)', line)
                date = match.group(1)
                # There is an annoying newline in the middle of the word version
                version = lines[m + 1].strip('ion').strip()
                break
            if line.startswith('Package: circuitikz'):
                match = re.search(r'circuitikz ([0-9/]+)', line)
                date = match.group(1)
                version = lines[m + 1].strip()
                break

        return date, version

    def find_circuitikz_version(self):

        content = r"""
        \documentclass[a4paper]{standalone}
        \usepackage{circuitikz}
        \begin{document}
        \end{document}"""

        tex_filename = tmpfilename('.tex')
        open(tex_filename, 'w').write(content)

        self.run(tex_filename)

        date, version = self.extract_circuitikz_version(tex_filename)

        self.cleanup(tex_filename)
        return date, version
