import os
import attr

@attr.s
class Paths(object):
    subdirs = attr.ib()
    root = attr.ib(default=os.getcwd())
    def __attrs_post_init__(self):
        for name in self.subdirs:
            path = os.path.join(self.root, name)
            object.__setattr__(self, name, path)

    def mkdirs(self):
        """Create `root` and `subdirs` if they don't already exist."""
        if not os.path.isdir(self.root):
            os.mkdir(self.root)
        for subdir in self.subdirs:
            path = self.__getattribute__(subdir)
            if not os.path.isdir(path):
                os.mkdir(path)


from .genome import Genome
from .species import Species
from .genbank import Genbank
from .metadata import Metadata
from .metadata import BioSample


# Figure out how to supress error output from this
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

__all__ = [Genome, Species, Genbank, Metadata, BioSample]
