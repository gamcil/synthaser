#!/usr/bin/env python3

import logging

from synthaser import settings
from synthaser.models import Serialiser

LOG = logging.getLogger(__name__)


def traverse_graph(graph, rules, domains, classifiers=None):
    """Traverses a rule graph and classifies a collection of domains.

    Each node is a dictionary with the schema:

        {
            "title": "Rule name",
            "children": [
                {
                    "title": Rule name",
                    "children": [ ... ],
                },
                ...
            ]
        }

    Rules are evaluated in order. If a rule is successfully evaluated,
    this function will recurse into any child rules, if any exist.

    Finally a classification list, containing the path of rules satisfied
    by the given domains, is returned.

    Args:
        graph (list, dict): Rule graph to traverse.
        rules (dict): Rule objects to evaluate on domains.
        domains (list): Domain objects to classify.
        classifiers (list): Current classifiers for a Domain collection.
    Returns:
        classifiers
    """
    if not classifiers:
        classifiers = []
    for node in graph:
        title = node["title"]
        rule = rules[title]
        if not rule.satisfied_by(domains):
            continue
        rule.rename_domains(domains)
        classifiers.append(title)
        children = node.get("children")
        if children:
            classifiers = traverse_graph(children, rules, domains, classifiers)
        return classifiers
    return classifiers


class RuleGraph(Serialiser):
    """A hierarchy of classification rules.

    The RuleGraph is used to classify synthases based on their domains.
    It stores Rule objects, as well as a directed graph controlling the
    order and hierarchy of classification.

    An example synthaser rule graph looks like this:

    ::

        [
            "Hybrid",
            {"PKS": ["HR-PKS", "PR-PKS", "NR-PKS"]},
            "NRPS"
        ]

    In this example, the "Hybrid" rule is evaluated first. If unsuccessful,
    the "PKS" rule is evaluated. If this is successful, synthaser recurses
    into child rules, in which case the "HR-PKS", "PR-PKS" and "NR-PKS" rules
    can be evaluated, and so on. Each rule name must have a corresponding
    entry in the rules attribute.

    Note that terminal leaves in the graph are placed in lists, whereas
    hierarchies are written as dictionaries of lists. This preserves rule
    order in Python, as well as preventing empty, unnecessary dictionaries
    at every level.

    Attributes:
        rules (dict): Collection of synthaser rules.
        graph (dict): Hierarchy of synthaser rules for classification.
    """

    def __init__(self, rules=None, graph=None):
        self.rules = rules if rules else {}
        self.graph = graph if graph else []

    @classmethod
    def from_dict(cls, d):
        return cls(
            rules={rule["name"]: Rule(**rule) for rule in d["rules"]},
            graph=d["hierarchy"]
        )

    def to_dict(self):
        return {
            "rules": [rule.to_dict() for rule in self.rules],
            "graph": self.graph
        }

    def classify(self, domains):
        return traverse_graph(self.graph, self.rules, domains)


class Rule:
    """A classification rule.

    Attributes:
        name (str): Name given to proteins satisfying this rule.
        domains (list): Domain types required to satisfy rule.
        filters (dict): Specific CDD families for each domain type.
        evaluator (str): Evaluatable rule satisfaction statement.
    """

    def __init__(
        self,
        name=None,
        renames=None,
        domains=None,
        filters=None,
        evaluator=None,
        **kwargs
    ):
        self.name = name if name else ""
        self.renames = renames if renames else []
        self.domains = domains if domains else []
        self.filters = filters if filters else []
        self.evaluator = evaluator if evaluator else ""

    def to_dict(self):
        return {
            "name": self.name,
            "renames": self.renames,
            "domains": self.domains,
            "filters": self.filters,
            "evaluator": self.evaluator
        }

    def evaluate(self, conditions):
        """Evaluates the rules evaluator string given evaluated conditions.

        Iterates backwards to avoid bad substitutions in larger (>=10) indices.
        e.g. "0 and 1 and ... and 13" --> "False and True and ... and True3"

        Arguments:
            conditions (list): Boolean values corresponding to domains in this rule.
        Returns:
            True if rule is satisfied, otherwise False.
        """
        evaluator = self.evaluator
        for idx, condition in reversed(list(enumerate(conditions))):
            evaluator = evaluator.replace(str(idx), str(condition))
        return eval(evaluator)

    def rename_domains(self, domains):
        """Renames domain types if substitutions are specified in the rule.

        The rename dictionary maps domain types to other domain types.
        For example, an ACP domain in a PKS matches the same PP-binding domain as
        a T domain in an NRPS, so to follow the naming convention the NRPS rule
        renames ACPs to Ts.

        Additionally, rename rules can be nested dicts to allow extra rules. For
        example, in a PKS-NRPS, the PP-binding domain in the NRPS module should
        be named T, not ACP. So, its rule is {'after': ['A', 'C'], 'to': 'T'};
        any ACP domains after the first A or C will be renamed T.
        """
        if not self.renames:
            return
        for rename in self.renames:
            flag = False
            for domain in domains:
                if not flag and domain.type in rename["after"]:
                    flag = True
                if flag and domain.type == rename["to"]:
                    domain.type = rename["to"]

    def valid_family(self, domain):
        """Checks a given domain matches a specified CDD family in the rule.

        If no families have been specified for the given domain type, this
        function will return True (i.e. any family of the type is accepted).

        This behaviour is controlled by the filters property of a synthaser rule.
        For example, to restrict a KS domain to certain CDD families:

        ::

            "filters": [
                "type": "KS",
                "domains": ["one", "two"]
            ]
        """
        for filt in self.filters:
            if filt["type"] == domain.type:
                return domain.accession in filt["domains"]
        return True

    def satisfied_by(self, domains):
        """Evaluates this rule against a collection of domains.

        Checks that:
        1) required domain types are represented in the supplied domains, and
        2) domains are of the desired CDD families, if any are specified.

        Placeholders in the evaluator string are then replaced by their
        respective booleans, and evaluated.

        Once a domain in the supplied domains has matched one in the rule, it
        cannot be matched to another in the rule. This enables rules based on
        counts of domains (e.g. multi-modular PKS w/ 2 KS domains).
        """
        LOG.debug("Evaluating %s against %s", self.name, [d.type for d in domains])
        seen = []
        conditions = []
        for rule_domain in self.domains:
            match = False
            for domain in domains:
                if domain in seen:
                    continue
                if domain.type == rule_domain and self.valid_family(domain):
                    seen.append(domain)
                    match = True
                    break
            conditions.append(match)
        return self.evaluate(conditions)


def classify(synthases, rule_file=None):
    """Classifies synthases based on defined rules.

    If no rule_file is provided, the packaged rules.json will be loaded by
    default.

    Arguments:
        synthases (list): Synthase objects to classify.
        rule_file (str): Path to custom classification rule file.
    """
    rule_file = rule_file or settings.RULE_FILE
    with open(rule_file) as fp:
        LOG.info("Loading rules: %s", fp.name)
        rg = RuleGraph.from_json(fp)
    for synthase in synthases:
        synthase.classification = rg.classify(synthase.domains)
