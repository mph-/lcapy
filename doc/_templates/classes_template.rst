{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
    :members:
    :show-inheritance:

    {% block methods %}

    {% if methods %}
    .. rubric:: {{ _('Methods') }}

    .. autosummary::
    {% for item in methods %}
        {%- if item not in inherited_members %}
            ~{{ name }}.{{ item }}
        {%- endif %}
    {%- endfor %}
    {% endif %}
    {% endblock %}

    {% block attributes %}
    {% if attributes %}
    .. rubric:: {{ _('Attributes') }}

    .. autosummary::
    {% for item in attributes %}
        {%- if item not in inherited_members %}
            ~{{ name }}.{{ item }}
        {%- endif %}
    {%- endfor %}
    {% endif %}
    {% endblock %}
