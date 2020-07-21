from collections import defaultdict


def group_synthases(synthases):
    """Group synthases by their classifications."""
    levels = defaultdict(list)
    for synthase in synthases:
        for level in synthase.classification:
            levels[level].append(synthase.header)
    return levels


def build_dict(path, d=None):
    """Recursively generates a dictionary of dictionaries from a list."""
    if not d:
        d = {}
    if path:
        key = path.pop(0)
        d[key] = build_dict(path, d)
    return d


def merge_dicts(a, b):
    """Recursively merges two dictionaries, allowing overlapping keys."""
    merged = dict(a)
    for key, value in b.items():
        if key in merged:
            merged[key] = merge_dicts(merged[key], value)
        else:
            merged[key] = value
    return merged


def get_classification_paths(synthases):
    """Determines the hierarchy of synthase classifications.

    This hierarchy is used when annotating the plot with classification bars.
    It should be used in conjunction with the per-classification synthase
    dictionary generated using group_synthases().
    """
    d = {}
    for synthase in synthases:
        classification = synthase.classification.copy()
        d = merge_dicts(d, build_dict(classification))
    return d


def iter_nested_keys(d, depth=0):
    """Iterates over all keys in a nested dictionary, reporting their depth.

    The depth indicates how deeply nested the yielded key is in the dictionary.
    It is used when annotating the plot to determine the position of the
    classification bars.
    """
    for key, values in d.items():
        yield key, depth
        if values:
            yield from iter_nested_keys(values, depth + 1)


def iter_annotation_groups(hierarchy):
    """Traverses hierarchy and iterates classification groups.

    Groups are reverse sorted by depth, such that annotations are drawn from more
    specific to less specific.
    """
    group = []
    for classification, depth in iter_nested_keys(hierarchy):
        if group and depth == 0:
            yield group
            group = []
        group.insert(0, {"classification": classification, "depth": depth})
    yield group


def get_annotation_groups(hierarchy):
    return [group for group in iter_annotation_groups(hierarchy)]
