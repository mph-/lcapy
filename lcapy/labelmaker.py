from .expr import Expr
from .latex import latex_format_label
from .engformatter import EngFormatter
import sympy as sym


class LabelMaker:

    def _format_expr(self, expr):

        return Expr(expr, cache=False).latex_math()

    def _format_value_units(self, value, units):

        return EngFormatter().latex_math(value, units)

    def _format_name(self, cpt_type, cpt_id):

        name = cpt_type
        subscript = cpt_id[1:] if cpt_id.startswith('_') else cpt_id

        if cpt_type == 'REL':
            name = r'\mathcal{R}'
        elif len(name) > 1:
            name = r'\mathrm{%s}' % name

        if subscript != '':
            if len(subscript) > 1:
                subscript = r'\mathrm{%s}' % subscript
            name = name + '_{%s}' % subscript
        return latex_format_label('$' + name + '$')

    def make(self, cpt, label_ports=False):

        # There are two possible labels for a component:
        # 1. Component name, e.g., R1
        # 2. Component value, expression, or symbol

        id_label = self._format_name(cpt.type, cpt.id)
        value_label = None

        if cpt.type == 'P' and not label_ports:
            id_label = None

        elif cpt.type in ('A', 'O', 'W') or id_label.find('#') != -1:
            id_label = None

        if cpt.type in ('A', 'S', 'SW', 'U'):
            value_label = ''

        unify = False
        if len(cpt.args):

            # TODO, extend for mechanical and acoustical components.
            units_map = {'V': 'V', 'I': 'A', 'R': '$\Omega$',
                         'C': 'F', 'L': 'H'}

            expr = cpt.args[0]
            if cpt.classname in ('Vstep', 'Istep'):
                expr = '(%s) * Heaviside(t)' % expr
                value_label = self._format_expr(expr)
            elif cpt.classname in ('Vs', 'Is'):
                value_label = self._format_expr(expr)
            elif cpt.classname == 'TF':
                expr = sym.sympify(expr)
                if expr.is_Pow and expr.args[1] == -1:
                    value_label = '%s:1' % (1 / expr)
                else:
                    value_label = '1:%s' % expr
            elif cpt.type in ('F', 'H') and len(cpt.args) > 1:
                # This is hard to give a reasonable label since the
                # control current is specified by a voltage source.
                # The user will have to override manually.
                expr = cpt.args[1]
                value_label = self._format_expr(expr)
            elif cpt.classname not in ('TP',):
                try:
                    # Handle things like 9/1000 that can occur
                    # when substituting cpt values.
                    value = float(sym.Rational(expr))
                    if cpt.type in units_map:
                        value_label = self._format_value_units(
                            value, units_map[cpt.type])
                    else:
                        value_label = self._format_expr(expr)

                except (ValueError, TypeError):
                    # This catches non numeric arg.
                    value_label = self._format_expr(expr)

            # Ensure labels are the same when the value is not specified.
            # This will prevent printing the name and value.
            unify = expr == cpt.type + cpt.id

            # Currently, we only annnotate the component with the
            # value, expression, or symbol.  If this is not specified,
            # it defaults to the component identifier.  Note, there
            # are some objects we do not want to label, such as wires
            # and ports.
        id_label = '' if id_label is None else latex_format_label(id_label)
        value_label = id_label if value_label is None \
            else latex_format_label(value_label)

        if unify:
            value_label = id_label

        return id_label, value_label
