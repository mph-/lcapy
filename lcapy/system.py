from os import system, path, remove, chdir, getcwd, stat
from subprocess import call, DEVNULL
import re
import platform

# System dependent functions

# Note, sometime in 2018 ImageMagick's convert program by default
# disallowed pdf file conversions.  A work-around is to edit
# /etc/ImageMagick-6/policy.xml to enable this conversion.

def tmpfilename(suffix=''):

    from tempfile import gettempdir, NamedTemporaryFile
    
    # Searches using TMPDIR, TEMP, TMP environment variables
    tempdir = gettempdir()
    
    filename = NamedTemporaryFile(suffix=suffix, dir=tempdir, 
                                  delete=False).name
    return filename


def convert_pdf_svg(pdf_filename, svg_filename):

    system('pdf2svg %s %s' % (pdf_filename, svg_filename))
    if not path.exists(svg_filename):
        raise RuntimeError('Could not generate %s with pdf2svg.  Is it installed?' % 
                           svg_filename)


def convert_pdf_png_convert(pdf_filename, png_filename, dpi=300):

    program = 'convert'    
    if platform.system() == 'Windows':
        program = 'magick convert'

    command = '%s -density %d %s %s' % (program, int(dpi), pdf_filename,
                                        png_filename)
    system(command)
        
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
            
    command = '%s -q -dQUIET -dSAFER -dBATCH -dNOPAUSE -dNOPROMPT -dMaxBitmap=500000000 -dAlignToPixels=0 -dGridFitTT=2 -sDEVICE=pngalpha -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -r%dx%d -sOutputFile=%s %s' % (program, int(dpi), int(dpi), png_filename, pdf_filename)

    system(command)

    
def convert_pdf_png_pdftoppm(pdf_filename, png_filename, dpi=300):

    system('pdftoppm -r %d  -png %s -thinlinemode shape > %s' %
           (int(dpi), pdf_filename, png_filename))
        
    if not path.exists(png_filename):
        raise RuntimeError('Could not generate %s with pdftoppm' % 
                           png_filename)    

    
def convert_pdf_png(pdf_filename, png_filename, dpi=300):

    try:
        convert_pdf_png_ghostscript(pdf_filename, png_filename, dpi)
    except:
        try:
            convert_pdf_png_convert(pdf_filename, png_filename, dpi)
        except:
            convert_pdf_png_pdftoppm(pdf_filename, png_filename, dpi)


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

    root, ext = path.splitext(tex_filename)
    dirname = path.dirname(tex_filename)
    baseroot = path.basename(root)
    cwd = getcwd()
    if dirname != '':
        chdir(path.abspath(dirname))
        
    call(['pdflatex', '-interaction', 'batchmode', '%s.tex' %
          baseroot], stderr=DEVNULL, stdout=DEVNULL)
    
    if dirname != '':
        chdir(cwd)            

    return root + '.pdf'


def run_dot(dotfilename, filename):

    base, ext = path.splitext(filename)
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

