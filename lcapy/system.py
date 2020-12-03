"""This module wraps system dependent programs for pdf generation, etc.

Copyright 2014--2020 Michael Hayes, UCECE

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


def checkexe(command):

    if isinstance(command, list):
        program = command[0]
        command = ' '.join(command)
    else:
        program = command.split(' ')[0]

    if False:
        print('Trying to run: %s' % command)
    if not hasexe(program):
        raise RuntimeError('%s is not installed' % program)


def run(command, stderr=DEVNULL, stdout=DEVNULL, shell=False):

    checkexe(command)
    call(command, stderr=stderr, stdout=stdout, shell=shell)
    
    
def tmpfilename(suffix=''):

    from tempfile import gettempdir, NamedTemporaryFile
    
    # Searches using TMPDIR, TEMP, TMP environment variables
    tempdir = gettempdir()
    
    filename = NamedTemporaryFile(suffix=suffix, dir=tempdir, 
                                  delete=False).name
    return filename


def convert_pdf_svg(pdf_filename, svg_filename):

    run(['pdf2svg', pdf_filename, svg_filename])
    if not path.exists(svg_filename):
        raise RuntimeError('Could not generate %s with pdf2svg.  Is it installed?' % 
                           svg_filename)


def convert_pdf_png_convert(pdf_filename, png_filename, dpi=300):

    program = 'convert'    
    if platform.system() == 'Windows':
        program = 'magick convert'

    run([program, '-density %d' % int(dpi), pdf_filename, png_filename])
        
    if not path.exists(png_filename):
        raise RuntimeError('Could not generate %s with convert' % 
                           png_filename)
    if stat(png_filename).st_size == 0:
        raise RuntimeError('Could not generate %s with convert, empty file' % 
                           png_filename)
    

def convert_pdf_png_ghostscript(pdf_filename, png_filename, dpi=300):

    program = 'gs'    
    if platform.system() == 'Windows':
        program = 'gswin32'
        if platform.machine().endswith('64'):
            program = 'gswin64'            

    run([program, '-q', '-dQUIET', '-dSAFER', '-dBATCH', '-dNOPAUSE',
         '-dNOPROMPT',  '-dMaxBitmap=500000000',  '-dAlignToPixels=0',
         '-dGridFitTT=2',  '-sDEVICE=pngalpha',  '-dTextAlphaBits=4',
         '-dGraphicsAlphaBits=4',  '-r%dx%d' % (int(dpi), int(dpi)),
         '-sOutputFile=' + png_filename, pdf_filename])

    
def convert_pdf_png_pdftoppm(pdf_filename, png_filename, dpi=300):

    root, ext = path.splitext(png_filename)
    args = ['pdftoppm', '-r %d' % int(dpi), '-png', '-thinlinemode shape', '-singlefile', pdf_filename, root]
    if False:
        # TODO, determine why this fails...
        run(args)
    else:
        from os import system
        checkexe(args)
        system(' '.join(args))
        
    if not path.exists(png_filename):
        raise RuntimeError('Could not generate %s with pdftoppm' % 
                           png_filename)    

    
def convert_pdf_png(pdf_filename, png_filename, dpi=300):

    conversions = (('ghostscript', convert_pdf_png_ghostscript),
                   ('convert', convert_pdf_png_convert),
                   ('pdftoppm', convert_pdf_png_pdftoppm))

    for conversion, func in conversions:
        try:
            return func(pdf_filename, png_filename, dpi)
        except:
            pass

    raise RuntimeError("""Could not convert pdf to png, tried: %s.  Check that one of these programs is installed.""" % ', '.join([conversion for conversion, func in conversions]))


def latex_cleanup(tex_filename, wanted_filename=''):

    root, ext = path.splitext(tex_filename)
    exts = ['.tex', '.aux', '.log', '.pdf']
    filenames = [root + ext for ext in exts]
    
    if wanted_filename in filenames:
        filenames.remove(wanted_filename)

    for filename in filenames:
        if path.exists(filename):
            remove(filename)


def run_latex(tex_filename):

    checkexe('pdflatex')
    
    root, ext = path.splitext(tex_filename)
    dirname = path.dirname(tex_filename)
    baseroot = path.basename(root)
    cwd = getcwd()
    if dirname != '':
        chdir(path.abspath(dirname))
        
    run(['pdflatex', '-interaction', 'batchmode', baseroot + '.tex'])
    
    if dirname != '':
        chdir(cwd)            

    return root + '.pdf'


def run_dot(dotfilename, filename):

    base, ext = path.splitext(filename)
    #run([which('dot'), '-T ' + ext[1:], '-o ' + filename, dotfilename])
    system('dot -T %s -o %s %s' % (ext[1:], filename, dotfilename))
    remove(dotfilename)            
    return


def circuitikz_version():

    content = r"""
    \documentclass[a4paper]{standalone}
    \usepackage{circuitikz}
    \begin{document}
    \end{document}"""

    tex_filename = tmpfilename('.tex')
    open(tex_filename, 'w').write(content)

    log_filename = tex_filename.replace('.tex', '.log')
    run_latex(tex_filename)

    with open(log_filename, 'rt') as logfile:
        lines = logfile.readlines()
        for m, line in enumerate(lines):
            if line.startswith('Package: circuitikz'):
                match = re.search(r'circuitikz ([0-9/]+)', line)
                date = match.group(1)
                version = lines[m + 1].strip()
                return date, version
    return None, None

