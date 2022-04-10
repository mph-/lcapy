from .quantity import Quantity


class TransferMixin(Quantity):

    quantity = 'transfer'
    quantity_label = 'Transfer function'
    quantity_units = ''
    is_always_causal = True
    # FIXME, non causal transfer functions can be constructed
    # but they won't result from circuit analysis.
    is_transfer = True
    is_ratio = True
