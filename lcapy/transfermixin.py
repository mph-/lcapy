from .quantity import Quantity

class TransferMixin(Quantity):

    quantity = 'transfer'
    quantity_label = 'Transfer function'
    units = ''
    is_always_causal = True
    is_transfer = True
