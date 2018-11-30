import attr
from pathlib import Path


@attr.s
class Paths(object):
    subdirs = attr.ib()
    root = attr.ib(default=Path.cwd(), validator=attr.validators.instance_of(Path))

    def __attrs_post_init__(self):
        for name in self.subdirs:
            path = self.root / name
            object.__setattr__(self, name, path)
        self.mkdirs()

    def mkdirs(self):
        """Create `root` and `subdirs` if they don't already exist."""
        self.root.mkdir(exist_ok=True)
        for subdir in self.subdirs:
            path = self.__getattribute__(subdir)
            path.mkdir(exist_ok=True)
