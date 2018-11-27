import os
from pathlib import Path

import attr


@attr.s
class Paths(object):
    subdirs = attr.ib()
    root = attr.ib(default=Path.cwd())

    def __attrs_post_init__(self):
        if not isinstance(self.root, Path):
            self.root = Path(self.root)
        for name in self.subdirs:
            path = self.root / name
            object.__setattr__(self, name, path)
        self.mkdirs()

    def mkdirs(self):
        """Create `root` and `subdirs` if they don't already exist."""
        if not os.path.isdir(self.root):
            os.mkdir(self.root)
        for subdir in self.subdirs:
            path = self.__getattribute__(subdir)
            if not os.path.isdir(path):
                os.mkdir(path)


from . import metadata as metadata
from .species import Species
from .genbank import Genbank
from .genome import Genome
from .metadata import BioSample
from .metadata import AssemblySummary

__all__ = [Genome, Species, Genbank, BioSample]


# Suppress error output with:
# wget -P ~/path/to/env/genbankqc/lib/ https://github.com/openwebos/qt/tree/master/lib/fonts
os.environ["QT_QPA_PLATFORM"] = "offscreen"
