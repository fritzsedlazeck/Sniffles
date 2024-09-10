from setuptools import setup, find_packages

setup(
    name='sniffles',
    version='2.4.1',
    packages=find_packages(),
    url='https://github.com/fritzsedlazeck/Sniffles',
    license='MIT',
    author='Moritz Smolka, Hermann Romanek',
    author_email='moritz.g.smolka@gmail.com, sniffles@romanek.at',
    description='A fast structural variation caller for long-read sequencing data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
