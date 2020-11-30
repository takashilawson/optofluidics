import setuptools
from os import path

dir = path.abspath(path.dirname(__file__))
with open(path.join(dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setuptools.setup(
    name='optofluidics',
    version='0.0.8',
    description='Optofluidics data processing',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/takashilawson/optofluidics',
    author='Takashi Lawson',
    author_email='takashilawson1@gmail.com',
    zip_safe=False,
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords='optofluidics spectroscopy hdf5 physics',
    packages=['optofluidics'],
    python_requires='>=3.5',
    install_requires=[
        'lmfit',
        'h5py',
        'pandas',
        'matplotlib',
        'scipy',
    ],
    include_package_data=True,
    )
