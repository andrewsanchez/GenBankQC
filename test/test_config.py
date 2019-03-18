from pathlib import Path
from tempfile import gettempdir

import pytest
from genbankqc import config


@pytest.fixture()
def simple_path():
    paths = config.Paths(root=gettempdir())
    return paths


def test_simple_path(simple_path):
    assert isinstance(simple_path.root, Path)
    assert simple_path.root.is_dir()


@pytest.fixture()
def nested():
    root = gettempdir()
    paths = config.Paths(root, subdirs=["a", "b", "c"])
    return paths


def test_nested(nested):
    assert isinstance(nested.root, Path)
    assert nested.root.is_dir()
    for d in nested.subdirs:
        assert getattr(nested, d).is_dir()


@pytest.fixture()
def nested_and_named():
    root = gettempdir()
    paths = config.Paths(root, subdirs=["a", "b", "c", ("special_name", "d")])
    return paths


def test_nested_and_named(nested_and_named):
    assert isinstance(nested_and_named.root, Path)
    assert nested_and_named.root.is_dir()
    for d in nested_and_named.subdirs:
        if isinstance(d, tuple):
            assert getattr(nested_and_named, d[0]).is_dir()
        else:
            assert getattr(nested_and_named, d).is_dir()
