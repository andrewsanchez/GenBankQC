import os

from . import config
from .metadata import Metadata
from .metadata import BioSample
from .metadata import AssemblySummary
from .species import Species
from .genbank import Genbank
from .genome import Genome

__all__ = ["Genome", "Species", "Genbank", "BioSample", "AssemblySummary", "Metadata"]


# Suppress error output with:
# wget -P ~/path/to/env/genbankqc/lib/fonts/ https://github.com/openwebos/qt/tree/master/lib/fonts
os.environ["QT_QPA_PLATFORM"] = "offscreen"
