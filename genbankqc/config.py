import attr
from pathlib import Path


@attr.s
class Paths(object):
    root = attr.ib(converter=Path)
    subdirs = attr.ib(default=[], validator=attr.validators.instance_of(list))

    def __attrs_post_init__(self):

        for subdir in self.subdirs:
            if isinstance(subdir, tuple):
                name, path = subdir
                path = self.root / path
                path.mkdir(exist_ok=True)
                object.__setattr__(self, name, path)
            else:
                path = self.root / subdir
                path.mkdir(exist_ok=True)
                object.__setattr__(self, subdir, path)
