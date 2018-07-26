from setuptools import setup, find_packages


setup(name='GenBankQC',
      version='0.2a',
      license='BSD-3-Clause',
      keywords='NCBI bioinformatics genomics',
      packages=find_packages(),
      python_requires='>=3.4',
      install_requires=[
          'click',
          'pandas',
          'pytest',
          'biopython',
          'retrying',
          'logbook',
      ],
      entry_points='''
      [console_scripts]
      genbankqc=genbankqc.__main__:cli
      ''',
      classifiers=[
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 3.6',
            'Operating System :: POSIX :: Linux',
            'License :: OSI Approved :: MIT License',
            'Intended Audience :: Science/Research',
            'Environment :: Console',
            'Development Status :: 3 - Alpha'])
