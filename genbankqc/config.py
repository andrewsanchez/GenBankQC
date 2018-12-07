import attr
from pathlib import Path


@attr.s
class Paths(object):
    root = attr.ib(converter=Path)
    subdirs = attr.ib(default=[], validator=attr.validators.instance_of(list))

    def __attrs_post_init__(self):
        self.root.mkdir(exist_ok=True)
        for subdir in self.subdirs:
            name = self.clean_path_name(subdir)
            path = self.root / subdir
            object.__setattr__(self, name, path)
        self.mkdirs()

    def mkdirs(self):
        """Create `root` and `subdirs` if they don't already exist."""
        for subdir in self.subdirs:
            name = self.clean_path_name(subdir)
            path = self.__getattribute__(name)
            path.mkdir(exist_ok=True)

    @staticmethod
    def clean_path_name(path):
        return path.strip(".")
