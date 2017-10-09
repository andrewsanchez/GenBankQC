from setuptools import find_packages, setup

package = 'genbank-qc'
version = '0.1'

setup(name=package,
      version=version,
      description="",
      packages=find_packages(),
      include_package_date=True,
      url='',
      install_requires=[
          'click',
          'numpy',
          'pandas',
          'scikit-bio',
          'biopython'
      ],
      entry_points='''
      [console_scripts]
      genbank-qc=genbank_qc.__main__:cli
      ''',
      )
