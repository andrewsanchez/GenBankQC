import os

from . import config
from . import metadata
from .species import Species
from .genbank import Genbank
from .genome import Genome
from .metadata import Metadata
from .metadata import BioSample
from .metadata import AssemblySummary

__all__ = ["Genome", "Species", "Genbank", "BioSample", "AssemblySummary", "Metadata"]


# Suppress error output with:
# wget -P ~/path/to/env/genbankqc/lib/fonts/ https://github.com/openwebos/qt/tree/master/lib/fonts
os.environ["QT_QPA_PLATFORM"] = "offscreen"
