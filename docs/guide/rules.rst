Creating custom rule sets
=========================

- synthaser originally created to analyse SM synthases
- designed to be completely modular, core functionality is agnostic to rules
- can therefore repurpose synthaser completely for any type of domains
  e.g. adding new domains, completely unrelated things

Domain rules
------------

- CD-search results in hundreds of overlapping domain hits
- visually, can see they segregate into distinct domain islands corresponding to a
  domain type
- synthaser programmatically determines these islands
- uses domain file to specify which domain hits are saved, what broader type they belong
  to etc
- define a handful of specific families per domain type, can then drastically lower hit
  quality thresholds without confusion (e.g. SAT)
- example KS domain rule
- how to get:
  - pssm length
  - bitscore threshold
- give to synthaser using --domain_file


Classification rules
--------------------

- Generally would like to classify sequences based on their domain architectures
- synthaser has a hierarchical classification system based on a rule file
- again, functionality is agnostic to actual chemical function, depends completely on
  your ruleset
- rules.json
  - examples
  - define one entry per rule
    - domains involved
    - name
    - specific domain families per type (filters)
    - evaluator
  - define hierarchy in a graph
    - written as arrays to preserve evaluation order
    - test
      - yes => if children, recurse and test
      - no => skip, go to next rule
- see synthaser.classify for example of creating rule graph file within python
