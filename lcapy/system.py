from os import system, path, remove, mkdir, chdir, getcwd
import re

# System dependent functions

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
        raise RuntimeError('Could not generate %s with pdf2svg' % 
                           svg_filename)


def convert_pdf_png(pdf_filename, png_filename, oversample=1):

    system('convert -density %d %s %s' %
           (oversample * 100, pdf_filename, png_filename))
    if path.exists(png_filename):
        return
        
    # Windows has a program called convert, try magick convert
    # for image magick convert.
    system('magick convert -density %d %s %s' %
           (oversample * 100, pdf_filename, png_filename))
    if path.exists(png_filename):
        return

    # Some Windows systems use im-convert.
    system('im-convert -density %d %s %s' %
           (oversample * 100, pdf_filename, png_filename))
    if path.exists(png_filename):
        return    

    raise RuntimeError('Could not generate %s with convert' % 
                       png_filename)


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
        
    system('pdflatex -interaction batchmode %s.tex' % baseroot)
    
    if dirname != '':
        chdir(cwd)            

    return root + '.pdf'


def run_dot(dotfilename, filename):

    base, ext = path.splitext(filename)
    system('dot -T %s -o %s %s' % (ext[1:], filename, dotfilename))
    remove(dotfilename)            
    return


def circuitikz_version():

    content = ('\\documentclass[a4paper]{standalone}\n'
               '\\usepackage{circuitikz}\n'
               '\\begin{document}\n'
               '\\end{document}\n')

    tex_filename = tmpfilename('.tex')
    open(tex_filename, 'w').write(content)

    log_filename = tex_filename.replace('.tex', '.log')
    run_latex(tex_filename)
    
    log = open(log_filename, 'r').read()
    latex_cleanup(tex_filename)    

    match = re.search(r'circuitikz ([0-9/]+)', log)
    if match is None:
        return None
    return match.group(1)
