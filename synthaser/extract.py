"""Extract domain/synthase sequences from synthaser results."""


import logging

from synthaser.models import SynthaseContainer


LOG = logging.getLogger(__name__)


def fasta(sequences):
    """Builds a FASTA record of extracted sequences.

    This function expects either a list of Synthase objects (mode='synthase')
    or a dictionary of extracted domains, keyed on synthase header (mode='domains').
    """
    # Entire Synthase objects were extracted
    if isinstance(sequences, list):
        return "\n".join(s.to_fasta() for s in sequences)

    # Domain sequences were extracted
    if isinstance(sequences, dict):
        return "\n".join(
            f">{header}_{index}\n{sequence}"
            for header, sequences in sequences.items()
            for index, sequence in enumerate(sequences)
        )


def write(sequences, prefix):
    """Writes extracted sequences to file(s) given some prefix."""
    for values in sequences.values():
        for k, synthases in values.items():
            with open(f"{prefix}{k}.faa", "w") as fp:
                LOG.info("  %s", fp.name)
                text = fasta(synthases)
                fp.write(text)


def extract(source, prefix, mode="domain", classes=None, types=None, families=None):
    """Extract domain or synthase sequences from a synthaser search session.

    Note that if mode='domains' and classes are specified, they will be used
    purely to filter the search session.
    """
    if not isinstance(source, SynthaseContainer):
        raise TypeError("Expected SynthaseContainer object")

    if mode == "domain":
        LOG.info("Extracting domain sequence(s) from synthase(s)")
        sequences = source.extract_domains(
            classes=classes,
            types=types,
            families=families,
            by="query",
        )
    elif mode == "synthase":
        LOG.info("Extracting synthase sequence(s)")
        sequences = source.extract_synthases(
            classes=classes,
            types=types,
            families=families,
        )
    else:
        raise ValueError("Expected 'domain' or 'synthase'")

    LOG.info("Writing sequences to file(s):")
    write(sequences, prefix)
