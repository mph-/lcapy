from .rcchecker import RcChecker
from .config import str_expr_map, latex_expr_map, pretty_expr_map, atoms_table
from .state import state

import sympy as sym


def imaginary_update(k, v):

    if v == 'i':
        str_expr_map[sym.I] = 'i'
        latex_expr_map[sym.I] = '\mathrm{i}'
        iunicode = '\u2148'
        atoms_table['ImaginaryUnit'] = iunicode
        pretty_expr_map[sym.I] = iunicode
    elif v == 'j':
        str_expr_map[sym.I] = 'j'
        latex_expr_map[sym.I] = '\mathrm{j}'
        junicode = '\u2149'
        atoms_table['ImaginaryUnit'] = junicode
        pretty_expr_map[sym.I] = junicode
    else:
        raise ValueError('symbols.imaginary must be i or j')


def heaviside_update(k, v):

    latex_expr_map[sym.Heaviside] = v
    pretty_expr_map[sym.Heaviside] = v


def sign_convention_update(k, v):

    state.sign_convention = v


def units_update(k, v):

    if k == 'units.loose':
        state.loose_units = v
    elif k == 'units.show':
        state.show_units = v
    elif k == 'units.abbreviate':
        state.abbreviate_units = v
    elif k == 'units.canonical':
        state.canonical_units = v
    elif k == 'units.check':
        state.check_units = v


def print_order_update(k, v):

    from .printing import printing_init
    printing_init(v)


c = RcChecker()

rcdefaults = {
    'sympy.solver' : ('DM', ('GJ', 'QR', 'CRAMER',
                             'GE', 'LU', 'ADJ', 'LDL', 'CH', 'DM')),
    'sympy.matrix.inverse' : ('DM', ('GE', 'LU', 'ADJ', 'LDL', 'CH', 'DM')),
    'sympy.print_order': ('lex', ('lex', 'grlex', 'grevlex'),
                          print_order_update),

    # Definition of H(0).  With H(0) = 0.5 then sgn(0) = 0 as expected
    # by SymPy and NumPy
    'functions.heaviside_zero' : (0.5, c.float),
    'functions.unitstep_zero' : (1, c.float),

    'symbols.imaginary' : ('j', ('i', 'j'), imaginary_update),
    'symbols.heaviside' : ('u', c.str, heaviside_update),

    'circuit.current_sign_convention': ('passive', ('passive', 'active'),
                                        sign_convention_update),

    'schematics.implicit_default' : ('sground', c.str),
    'schematics.autoground_default' : ('sground', c.str),

    'schematics.draw_nodes': ('primary', ('all', 'none', 'primary', 'connections')),
    'schematics.label_nodes': ('primary', ('all', 'alpha', 'none', 'primary')),
    'schematics.label_values': (True, c.bool),
    'schematics.label_ids': (True, c.bool),
    'schematics.label_style': ('aligned',
                               ('aligned', 'stacked', 'split', 'value', 'name')),
    'schematics.label_value_style': ('eng3', c.str),
    'schematics.label_flip': (False, c.bool),
    'schematics.annotate_values': (False, c.bool),
    'schematics.anchor': ('south east', c.str),
    'schematics.autoground': ('none', c.str),
    'schematics.scale': (1.0, c.float),
    'schematics.dpi': (300, c.int),
    'schematics.cpt_size': (1.5, c.float),
    'schematics.node_spacing': (2.0, c.float),
    'schematics.help_lines': (0.0, c.float),
    'schematics.style': ('american', ('american', 'british', 'european')),
    'schematics.voltage_dir': ('RP', ('RP', 'EF')),

    'os.linux.ghostscript' : ('gs', c.str),
    'os.macos.ghostscript' : ('gs', c.str),
    'os.windows32.ghostscript' : ('gswin32c', c.str),
    'os.windows64.ghostscript' : ('gswin64c', c.str),

    'os.linux.convert' : ('convert', c.str),
    'os.macos.convert' : ('convert', c.str),
    'os.windows32.convert' : ('magick convert', c.str),
    'os.windows64.convert' : ('magick convert', c.str),

    'os.linux.pdf2svg' : ('pdf2svg', c.str),
    'os.macos.pdf2svg' : ('pdf2svg', c.str),
    'os.windows32.pdf2svg' : ('pdf2svg', c.str),
    'os.windows64.pdf2svg' : ('pdf2svg', c.str),

    'os.linux.pdflatex' : ('pdflatex', c.str),
    'os.macos.pdflatex' : ('pdflatex', c.str),
    'os.windows32.pdflatex' : ('pdflatex', c.str),
    'os.windows64.pdflatex' : ('pdflatex', c.str),

    'os.linux.dot' : ('dot', c.str),
    'os.macos.dot' : ('dot', c.str),
    'os.windows32.dot' : ('dot', c.str),
    'os.windows64.dot' : ('dot', c.str),

    # Allow 1 + s, etc.
    'units.loose': (True, c.bool, units_update),
    # Print units with expression.
    'units.show': (False, c.bool, units_update),
    # Print units in abbreviated form.
    'units.abbreviate': (True, c.bool, units_update),
    # Print units in canonical form.
    'units.canonical': (False, c.bool, units_update),
    # This cannot be enabled without loose units (e.g.,
    # s + 1, or s * 'R' * 'C' + 1 would fail).
    'units.check': (True, c.bool, units_update),
}
