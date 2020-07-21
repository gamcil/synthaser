.. _models_module:

:mod:`synthaser.models`
-------------------------

This module stores the classes used throughout synthaser.

The Domain class represents a conserved domain hit. It stores the broader domain type,
the specific conserved domain profile name (from CDD), as well as its position in its
parent synthase sequence and score from the search. It also provides methods for slicing
the corresponding sequence and serialisation. We can instantiate a `Domain` object like
so:

>>> from synthaser.models import Domain
>>> domain = Domain(
...     type='KS',
...     domain='PKS_KS',
...     start=756,
...     end=1178,
...     evalue=0.0,
...     bitscore=300
... )

and get its sequence given the parent `Synthase` object sequence:

>>> domain.slice(synthase.sequence)
'MPIAVGM..'

Likewise, the `Synthase` class stores information about a synthase, including its name,
amino acid sequence, `Domain` instances and its classification. It also
contains methods for generating the domain architecture, extraction of domain sequences
and more. For example, we can instantiate a new `Synthase` object like so:

>>> from synthaser.models import Synthase
>>> synthase = Synthase(
...     header='SEQ001.1',
...     sequence='MASGTC...',
...     domains=[
...         Domain(type='KS'),
...         Domain(type='AT'),
...         Domain(type='DH'),
...         Domain(type='ER'),
...         Domain(type='KR'),
...         Domain(type='ACP'),
...     ],
... )

Then, we can generate the domain architecture:

>>> synthase.architecture
'KS-AT-DH-ER-KR-ACP'

Or extract all of the domain sequences:

>>> synthase.extract_domains()
{
    "KS_0": "MPIAVGM...",
    "AT_0": "VFTGQGA...",
    "DH_0": "DLLGVPV...",
    "ER_0": "DVEIQVS...",
    "KR_0": "IAENMCS...",
    "ACP_0": "ASTTVAQ..."
}

The object can also be serialised to JSON (note the Domain object works the same way):

>>> js = synthase.to_json()
>>> with open('synthase.json', 'w') as handle:
...     handle.write(js)

and subsequently loaded from JSON:

>>> with open('synthase.json') as handle:
...     synthase = Synthase.from_json(handle)

This will internally convert the Synthase object, as well as any Domain objects it
contains, to dictionaries, before converting to JSON using the builtin json library
and writing to file. When loading up from JSON, this process is reversed, and the
entries in the file are converted back to Python objects.

.. automodule:: synthaser.models
        :members:
