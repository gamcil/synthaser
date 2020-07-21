.. _classify_module:

:mod:`synthaser.classify`
-------------------------

This module contains the logic for classifying synthase objects based on
user-defined rules.

To classify a collection of sequences, use the classify function:

>>> from synthaser.classify import classify
>>> classify(my_sequences)

A custom classification rule file can be provided to this function like so:

>>> classify(my_sequences, rule_file="my_rules.json")

Briefly, rule files should contain:

1) Rule entries, specifying the domain combinations required to satisfy them
2) Rule graph, encoding the hierarchy and order in which rules are evaluated

Alternatively, you could build a `RuleGraph` object in Python, e.g.:

>>> from synthaser.classify import Rule, RuleGraph
>>> one = Rule(name="Rule 1", domains=["D1", "D2"], evaluator="0 and 1")
>>> two = Rule(name="Rule 2", domains=["D3", "D4", "D5"], evaluator="(0 and 1) or 2")
>>> three = Rule(name="Rule 3", domains=["D6", "D7"], evaluator="0 or 1")
>>> graph = [
...     "Rule 1",
...     {
...         "Rule 2": ["Rule 3"]
...     }
... ]
>>> rg = RuleGraph(rules=[one, two, three], graph=graph)

And then save it to a file:

>>> with open("my_rules.json", "w") as fp:
...     rg.to_json(fp)

This `RuleGraph` object can directly classify `Synthase` objects:

>>> rg.classify(my_sequences)

For further explanation of rule files, refer to the documentation.

.. automodule:: synthaser.classify
        :members:
