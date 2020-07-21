.. _plot_module:

:mod:`synthaser.plot`
-------------------------

This modules handles the construction of `synthaser` visualisations.

`synthaser` plots are implemented as HTML documents with embedded visualisations
created using the D3 JavaScript library. Files used in this process are all stored separately
inside the `plot` folder in the source code (i.e. separate CSS, JavaScript, HTML).

This module provides two ways to construct a `synthaser` plot from a collection of
`Synthase` objects: hosting using Python's built in `socketserver` module, or generating
a completely static HTML file containing all code elements necessary for the plot.

1. Get necessary data

>>> data = plot.get_data(synthases)

This generates a dictionary which can be converted easily to JSON and given to the
JavaScript visualisation. It contains a dictionary of the `Synthase` objects, the order in
which they should be drawn, a dictionary of `Synthase` objects grouped by their
classifications, and the annotation group hierarchy generated using
`grouping.get_annotation_groups()`.

2. a) Dynamically serve plots using `socketserver`

>>> plot.serve_html(data)

This will host the visualisation on some randomly chosen open port on localhost. It
requires a keyboard interrupt to stop serving the plot.

2. b) Generate static HTML file

>>> plot.save_html(data, "myoutput.html")

This will generate a completely static HTML file containing all code elements required
to render the plot (JavaScript and CSS embedded directly in the HTML file). This can
then be moved/copied anywhere you would like.

This workflow is mirrored by the `plot_synthases` function:

>>> plot.plot_synthases(synthases, output="myoutput.html")


.. automodule:: synthaser.plot
        :members:
