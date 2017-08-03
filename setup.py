from setuptools import setup, find_packages

package = 'genbankfilter'
version = '0.1'

setup(name=package,
      version=version,
      description="",
      packages=find_packages(),
      include_package_date=True,
      url='',
      install_requires=[
          'click',
          'ncbitk',
          'biopython',
          'matplotlib'
      ],
      entry_points='''
      [console_scripts]
      gbfilter=genbankfilter.__main__:cli
      ''',
      )
