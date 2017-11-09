from setuptools import find_packages, setup


setup(name='genbank-qc',
      version='0.1a3',
      license="MIT",
      url="https://github.com/andrewsanchez/genbank-qc",
      description="Automated quality control for Genbank genomes.",
      author='Andrew Sanchez',
      author_email='inbox.asanchez@gmail.com',
      keywords='NCBI bioinformatics',
      packages=find_packages(),
      python_requires='>=3.4',
      install_requires=[
          'click',
          'pandas',
          'pytest',
          'biopython'],
      entry_points='''
      [console_scripts]
      genbank-qc=genbank_qc.__main__:cli
      ''',
      classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 3.6',
            'Operating System :: POSIX :: Linux',
            'License :: OSI Approved :: MIT License',
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            'Development Status :: 3 - Alpha'])
