from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pyticle_tracker',
      version='0.5',
      description='Lagrangian particle tracking for FVCOM',
      long_description=readme(),
      url='https://github.com/moflaher/pyticle_tracker',
      author='Mitchell O\'Flaherty-Sproul, Kody Crowell',
      author_email='073208o@acadiau.ca, 119865c@acadiau.ca',
      packages=find_packages(),
      zip_safe=False)
