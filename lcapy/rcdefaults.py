from .rcchecker import RcChecker
from .config import str_expr_map, latex_expr_map, pretty_expr_map, atoms_table
from .state import state

import sympy as sym


def imaginary_update(v):

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


def heaviside_update(v):

    latex_expr_map[sym.Heaviside] = v
    pretty_expr_map[sym.Heaviside] = v


def sign_convention_update(v):

    state.sign_convention = v


c = RcChecker()

rcdefaults = {
    'sympy.solver' : ('DM', c.str),
    'sympy.matrix.inverse' : ('DM', c.str),
    'sympy.matrix.inverse-fallback' : ('ADJ', c.str),

    'functions.heaviside-zero' : (0.5, c.float),
    'functions.unitstep-zero' : (0, c.float),

    'schematics.symbols.implicit' : ('sground', c.str),
    'schematics.symbols.autoground' : ('sground', c.str),

    'symbols.imaginary' : ('j', ('i', 'j'), imaginary_update),
    'symbols.heaviside' : ('u', c.str, heaviside_update),

    'current.sign-convention': ('passive', ('passive', 'active'),
                                sign_convention_update),

    'schematics.draw_nodes': ('primary', ('all', 'none', 'primary', 'connections')),
    'schematics.label_nodes': ('primary', ('all', 'alpha', 'none', 'primary')),
    'schematics.label_values': (True, c.bool),
    'schematics.label_ids': (True, c.bool),
    'schematics.label_style': ('aligned', ('aligned', 'stacked', 'split', 'value', 'name')),
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
    'schematics.voltage_dir': ('RP', ('RP', 'EF'))
}
