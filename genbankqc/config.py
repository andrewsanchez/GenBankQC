import attr
from pathlib import Path


@attr.s
class Paths(object):
    subdirs = attr.ib()
    root = attr.ib(default=Path.cwd())

    def __attrs_post_init__(self):
        if not isinstance(self.root, Path):
            self.root = Path(self.root)
        for subdir in self.subdirs:
            name = self.clean_path_name(subdir)
            path = self.root / subdir
            object.__setattr__(self, name, path)
        self.mkdirs()

    def mkdirs(self):
        """Create `root` and `subdirs` if they don't already exist."""
        self.root.mkdir(exist_ok=True)
        for subdir in self.subdirs:
            name = self.clean_path_name(subdir)
            path = self.__getattribute__(name)
            path.mkdir(exist_ok=True)

    @staticmethod
    def clean_path_name(path):
        return path.strip('.')
