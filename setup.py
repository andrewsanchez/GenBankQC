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
          'pandas',
          'biopython',
          'matplotlib',
          'scipy',
          'ncbitk',
      ],
      entry_points='''
      [console_scripts]
      gbf=genbankfilter.__main__:cli
      ''',
      )
