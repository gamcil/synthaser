.. _grouping_module:

:mod:`synthaser.grouping`
-------------------------

This module contains some functions used for grouping `Synthase` objects
by their classifications.

This is used primarily when grouping sequences for the purpose of annotation in the plot
(i.e. grouping Synthases of like classification, at each level in the classification
hierarchy). Since annotations need to be drawn from more specific to less specific,
this module generates groups in reverse.

Given a collection of classified `Synthase` objects, a basic workflow using this
module might be:

1. Build a dictionary of synthase headers grouped by classification:

>>> levels = group_synthases(synthases)
>>> levels
defaultdict(<class 'list'>, {'PKS': ['seq1', 'seq2', ...], 'HR-PKS': ['seq1', ...]})

2. Determine the hierarchy of synthase classifications in your synthases.

>>> hierarchy = get_classification_paths(synthases)
>>> hierarchy
{'PKS': {'Type I': {'Non-reducing': {}, 'Highly-reducing': {}, 'Partially-reducing': {}}}, 'Hybrid': {}}

Note, this is agnostic to our rule files - the rule hierarchy here is built solely
from what is stored in each `Synthase` object. This also means there should be no
redundant classifications.

3. Build an array of annotation groups, each in drawing (i.e. reverse) order.

>>> groups = get_annotation_groups(hierarchy)
>>> groups
[
  [
    {'classification': 'Partially-reducing', 'depth': 2},
    {'classification': 'Highly-reducing', 'depth': 2},
    {'classification': 'Non-reducing', 'depth': 2},
    {'classification': 'Type I', 'depth': 1},
    {'classification': 'PKS', 'depth': 0}
  ],
  [{'classification': 'Hybrid', 'depth': 0}],
]

Since annotations are drawn from more to less specific, and each classification is drawn
at some offset to the previous one, we need some way of differentiating their level -
hence the `depth` property.

.. automodule:: synthaser.grouping
        :members:
