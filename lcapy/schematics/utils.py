def check_boolean(value):

    if value not in (True, False, 'none', None, 'true', 'false'):
        raise ValueError('Unexpected Boolean value %s' % value)
    return value in (True, 'true')
